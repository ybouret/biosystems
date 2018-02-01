#include <Servo.h>
#include <medium.h>

////////////////////////////////////////////////////////////////////////////////
//
// Parameters for the control of the experiment
//
////////////////////////////////////////////////////////////////////////////////

//! PROBE_DELAY in seconds
/**
  the global sensing delay
*/
#define PROBE_DELAY       1.0f


//! PROBE_FREQUENCY in Hz
#define PROBE_FREQUENCY   40.0f

//! OUTPUT_FREQUENCY in Hz
#define OUTPUT_FREQUENCY  5.0f

////////////////////////////////////////////////////////////////////////////////
//
// Board setup
//
////////////////////////////////////////////////////////////////////////////////

#define servo_pin 9
#define pinFn     1
#define pinFt     4
#define baudRate  115200

////////////////////////////////////////////////////////////////////////////////
//
// Medium management: I/O
//
////////////////////////////////////////////////////////////////////////////////
static Medium medium;

////////////////////////////////////////////////////////////////////////////////
//
// Fish possible states
//
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
//
// States names
//______________________________________________________________________________
enum State
{
  PRIMING,        //!< accumulate data in memory
  RUNNING,        //!< start acting according to model
  OUT_OF_MEMORY   //!< parameters are too demanding!
};

//______________________________________________________________________________
//
// default state
//______________________________________________________________________________
static State status = PRIMING;

////////////////////////////////////////////////////////////////////////////////
//
// Servo management
//
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
//
// Servo parameters
//______________________________________________________________________________
static Servo servo;
static float         servo_angle_mini = 0.0f;
static float         servo_angle_maxi = 180.0f;
static float         servo_angle_init = 90.0f;
static unsigned long servo_pwm_lower  = 750;
static unsigned long servo_pwm_upper  = 2250;
static float         servo_last_time  = 0.0f;
static float         servo_period     = 10.0f;

//______________________________________________________________________________
//
// Servo_attach(), to be called in setup()
//______________________________________________________________________________
static void Servo_attach()
{
  servo.detach();
  servo.attach(servo_pin, servo_pwm_lower, servo_pwm_upper);
  servo_last_time = Medium_GetCurrentTime();
}

////////////////////////////////////////////////////////////////////////////////
//
// Forces
//
////////////////////////////////////////////////////////////////////////////////
static float FRange  = 5.0f;  //!< force range in Volts
static float FnInit  = 2.50f; //!< default normal     force reading (means 0 force)
static float FtInit  = 2.50f; //!< default tangential force reading (means 0 force)

#define READ_FORCE_ON(PIN) ( (FRange * ( (float) analogRead(PIN) ))/1023.0f )
#define READ_FORCE_N() READ_FORCE_ON(pinFn)
#define READ_FORCE_T() READ_FORCE_ON(pinFt)


////////////////////////////////////////////////////////////////////////////////
//
// Fish memory management
//
////////////////////////////////////////////////////////////////////////////////

//______________________________________________________________________________
//
// fields to be kept in memory
//______________________________________________________________________________

#define MEMORY_FIELDS 3

#define MEMORY_ANGLE  0
#define MEMORY_FN     1
#define MEMORY_FT     2

//______________________________________________________________________________
//
// memory timings management
//______________________________________________________________________________
static float               memory_delta_t     = 0.0f;     //!< mininal time between nodes
static const float PROGMEM memory_extra_delay = 10.0e-3f; //!< extra time in seconds
static float               memory_last_time   = 0.0f;     //!< time control

//______________________________________________________________________________
//
// a memory node, to store data into a list
//______________________________________________________________________________
struct memoryNode
{
  memoryNode *next;                  //!< for list
  memoryNode *prev;                  //!< for list
  float       time;                  //!< recorded time
  float       fields[MEMORY_FIELDS]; //!< space for fields
  int         quality;               //!< quality: do we respect delta_t

  //! probe current time and fields
  /**
     \param old if not NULL, compute quality
  */
  inline void probe(const memoryNode *old) throw()
  {
    time = Medium_GetCurrentTime();
    fields[MEMORY_ANGLE] = servo.read();
    fields[MEMORY_FN]    = READ_FORCE_N();
    fields[MEMORY_FT]    = READ_FORCE_T();
    if (NULL != old)
    {
      const float real_dt = time - old->time;
      quality = (int)floorf( (100.0f * memory_delta_t / real_dt) + 0.5f );
    }
    else
    {
      quality = 0;
    }

  }

  //! code factorisation to copy fields+quality
  inline void copy_data_from( const memoryNode *node ) throw()
  {
    memcpy(fields, node->fields, sizeof(fields));
    quality = node->quality;
  }
};

typedef Medium::ListOf<memoryNode> _memoryList; //!< base class for memory
typedef Medium::PoolOf<memoryNode>  memoryPool; //!< to handle nodes

//______________________________________________________________________________
//
//! doubly linked list of memory nodes
//______________________________________________________________________________
class memoryList : public _memoryList
{
  public:
    //! default constructor, initialize
    inline explicit memoryList() throw() :
      _memoryList(),
      cache(),
      nodes(NULL)
    {

    }

    //! destructor, free memory
    inline virtual ~memoryList() throw()
    {
      if (nodes) free(nodes);
    }

    //! allocate nodes and put them in the cache
    void setup(const unsigned long num_nodes) throw()
    {
      // allocated memory
      const unsigned long required = num_nodes * sizeof(memoryNode);
      nodes = (memoryNode *)malloc(required);
      if (!nodes)
      {
        status = OUT_OF_MEMORY; // will do nothing but complain!
        return;
      }

      // initialize cache
      memset(nodes, 0, required);
      for (unsigned long i = 0; i < num_nodes; ++i)
      {
        cache.store( nodes + i );
      }
    }

    //! populate the memory while some nodes are in the cache
    inline void prime() throw()
    {
      if (cache.size > 0)
      {
        memoryNode *node = cache.query();
        node->probe(head);
        push_front(node);
        node->fields[MEMORY_ANGLE] = 90.0f + 30.0f * Medium::SineWave(node->time, PROBE_DELAY / 2.0f);
      }
      else
      {
        status = RUNNING;
      }
    }

    //! remove the last node and make a new first node from it
    inline void update()
    {
      memoryNode *node = pop_back();
      node->probe(head);
      push_front(node);
#if 0
      Serial.print(F("t="));  Serial.print(node->time);
      Serial.print(F(" Q=")); Serial.print(node->quality);
      Serial.println(F(""));
#endif
    }

    //! interpolate a node from the memory@t-PROBE_DELAY
    inline void interpolate(const float t, memoryNode *node ) const
    {
      const float time = (node->time = t - PROBE_DELAY);

      const memoryNode *curr = tail;
      if (time <= curr->time)
      {
        // shouldn't happen, the last node is too young!
        node->copy_data_from(curr);
      }
      else
      {
        // bracket the target time
        const memoryNode *prev = curr->prev;
        while (NULL != prev)
        {
          const float c_time = curr->time;
          const float p_time = prev->time;
          if ( (c_time <= time) && (time <= p_time) )
          {
            // linear interpolation of all fields
            const float  fac = (time - c_time) / (p_time - c_time);
            const float *fc  = curr->fields;
            const float *fp  = prev->fields;
            float       *fn  = node->fields;
            for (unsigned i = 0; i < MEMORY_FIELDS; ++i)
            {
              const float yc = fc[i];
              const float yp = fp[i];
              fn[i] = yc + fac * (yp - yc);
            }
            node->quality = (curr->quality + prev->quality) / 2;
            return;
          }
          curr = prev;
          prev = curr->prev;
        }
        // shouldn't happen, the first node is too old!!
        node->copy_data_from(curr);
      }
    }

  private:
    memoryList(const memoryList &);
    memoryList&operator=(const memoryList &);

    memoryPool  cache; //!< the nodes cache
    memoryNode *nodes; //!< memory to hold the nodes

};

//______________________________________________________________________________
//
//! the memory...
//______________________________________________________________________________
static memoryList memory;

//______________________________________________________________________________
//
//! to be called in setup()
//______________________________________________________________________________
static void Memory_setup()
{
  // compute #nodes
  unsigned long num_nodes = 1 + (unsigned long)( ceilf(PROBE_DELAY * PROBE_FREQUENCY) );
  if (num_nodes < 2)
  {
    num_nodes = 2;
  }

  // prepare history
  Serial.print(F("Allocating #nodes=")); Serial.println(num_nodes);
  memory.setup(num_nodes);

  // compute effective delay
  const float total_delay = PROBE_DELAY + memory_extra_delay;
  memory_delta_t = total_delay / (num_nodes - 1);
  Serial.print(F("memory_delta_t=")); Serial.println(memory_delta_t, 6);

  // prepare loop timing
  memory_last_time = Medium_GetCurrentTime();
}

//______________________________________________________________________________
//
//! memory priming function
//______________________________________________________________________________
static void Memory_priming()
{
  const float t = Medium_GetCurrentTime();
  if (t - memory_last_time >= memory_delta_t)
  {
    Serial.println(F("...priming"));
    memory.prime();
    memory_last_time = t;
  }
}

//______________________________________________________________________________
//
//! memory running function
//______________________________________________________________________________
static void Memory_loop()
{
  const float t = Medium_GetCurrentTime();
  if (t - memory_last_time >= memory_delta_t)
  {
    memory.update();
    memory_last_time = t;
  }
}

//______________________________________________________________________________
//
//! servo control for main loop
//______________________________________________________________________________
static void Servo_loop()
{
  const float t = Medium_GetCurrentTime();
  memoryNode  node;
  memory.interpolate(t, &node);
  //servo.write( 180.0f * Medium::Triangle(t, servo_period) );
  const float angle = node.fields[MEMORY_ANGLE];
  servo.write(angle);
  if (t - servo_last_time >= 1.0 / OUTPUT_FREQUENCY)
  {
    Serial.print("t=");      Serial.print(t);
    Serial.print(" angle="); Serial.print(angle);
    Serial.print(" Q=");     Serial.print(node.quality);
    Serial.println("");
    servo_last_time = t;
  }
}

////////////////////////////////////////////////////////////////////////////////
//
// global setup
//
////////////////////////////////////////////////////////////////////////////////
void setup()
{
  // I/O
  Serial.begin(baudRate);
  while (!Serial);

  // all the setup functions
  Servo_attach();
  Memory_setup();
}

////////////////////////////////////////////////////////////////////////////////
//
// input processing
//
////////////////////////////////////////////////////////////////////////////////

//! table of commands into flash memory
static  const char * const PROGMEM commands[]  =
{
  "period"
};

//! command aliases
#define __PERIOD commands[0]

static void processInput()
{
  const unsigned numWords = medium.splitInput();
  if (2 == numWords)
  {
    const char *cmd = medium[0];
    if ( Medium_streq(__PERIOD, cmd) )
    {
      const float tmp = atof(medium[1]);
      if (tmp > 0)
      {
        servo_period = tmp;
      }
    }
  }

  goto END_INPUT;
END_INPUT:
  medium.resetInput();
}

////////////////////////////////////////////////////////////////////////////////
//
// global loop
//
////////////////////////////////////////////////////////////////////////////////
void loop()
{
  const float t = Medium_GetCurrentTime();
  switch (status)
  {
    case PRIMING:
      Memory_priming();
      break;

    case RUNNING:
      Memory_loop();
      Servo_loop();
      break;

    case OUT_OF_MEMORY:
      if (t - servo_last_time >= 1)
      {
        Serial.println(F("#Out of memory..."));
        servo.write(servo_angle_init);
        servo_last_time = t;
      }
      break;
  }

  // I/O
  if (medium.inputCompleted() )
  {
    processInput();
  }
}

////////////////////////////////////////////////////////////////////////////////
//
// global I/O
//
////////////////////////////////////////////////////////////////////////////////
void serialEvent()
{
  medium.serialEventCallback();
}


