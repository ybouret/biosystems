#include <Servo.h>
#include <medium.h>

//! PROBE_DELAY in seconds
#define PROBE_DELAY       1.0f

//! PROBE_FREQUENCY in Hz
#define PROBE_FREQUENCY   5.0f

// Board setup
#define servo_pin 9
#define pinFn     1
#define pinFt     4
#define baudRate  115200

// Medium management
static Medium medium;

enum State
{
  PRIMING,
  RUNNING,
  OUT_OF_MEMORY
};

static State status = PRIMING;

// Servo management
static Servo servo;
static float servo_angle_mini = 0.0f;
static float servo_angle_maxi = 180.0f;
static float servo_angle_init = 90.0f;
static unsigned long servo_pwm_lower = 750;
static unsigned long servo_pwm_upper = 2250;
static float servo_last_time = 0.0f;

static void Servo_attach()
{
  servo.detach();
  servo.attach(servo_pin, servo_pwm_lower, servo_pwm_upper);
  servo_last_time = Medium_GetCurrentTime();
}

// Forces
static float FRange  = 5.0f;  //!< force range in Volts
static float FnInit  = 2.50f; //!< initial normal force reading (means 0 force)
static float FtInit  = 2.50f; //!< initial tangential force reading (means 0 force)

#define READ_FORCE_ON(PIN) ( (FRange * ( (float) analogRead(PIN) ))/1023.0f )
#define READ_FORCE_N() READ_FORCE_ON(pinFn)
#define READ_FORCE_T() READ_FORCE_ON(pinFt)


// fish memory setup
#define MEMORY_FIELDS 3
#define MEMORY_ANGLE  0
#define MEMORY_FN     1
#define MEMORY_FT     2


static float               memory_delta_t     = 0.0f;
static const float PROGMEM memory_extra_delay = 10.0e-3f; //!< extra time in seconds
static float               memory_last_time   = 0.0f;

struct memoryNode
{
  memoryNode *next;
  memoryNode *prev;
  float       time;
  float       fields[MEMORY_FIELDS];
  int         quality;
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
  }
};

typedef Medium::ListOf<memoryNode> _memoryList;
typedef Medium::PoolOf<memoryNode>  memoryPool;

class memoryList : public _memoryList
{
  public:
    inline explicit memoryList() throw() :
      _memoryList(),
      cache(),
      nodes(NULL)
    {

    }

    inline virtual ~memoryList() throw()
    {
      if (nodes) free(nodes);
    }

    // allocate nodes and put them in the cache
    void setup(const unsigned long num_nodes)
    {
      const unsigned long required = num_nodes * sizeof(memoryNode);
      nodes = (memoryNode *)malloc(required);
      if (!nodes)
      {
        status = OUT_OF_MEMORY;
        return;
      }
      Serial.println(F("...allocated"));
      memset(nodes, 0, required);
      for (unsigned long i = 0; i < num_nodes; ++i)
      {
        cache.store( nodes + i );
      }
    }

    inline void prime()
    {
      if (cache.size > 0)
      {
        memoryNode *node = cache.query();
        node->probe(head);
        push_front(node);
      }
      else
      {
        status = RUNNING;
      }
    }

    inline void update()
    {
      memoryNode *node = pop_back();
      node->probe(head);
      push_front(node);
      Serial.print(F("t="));  Serial.print(node->time);
      Serial.print(F(" Q=")); Serial.print(node->quality);
      Serial.println(F(""));
    }

    inline void interpolate(const float t, memoryNode *node ) const
    {
      const float time = (node->time = t - PROBE_DELAY);
      //Serial.print(head->time); Serial.print(" "); Serial.print(time); Serial.print(" "); Serial.println(tail->time);

      const memoryNode *curr = tail;
      if (time <= curr->time)
      {
        // shouldn't happen, the last node is too young!
        memcpy(node->fields, curr->fields, sizeof(memoryNode::fields));
        node->quality = curr->quality;
        //Serial.println(">");
      }
      else
      {
        // bracket the target time
        const memoryNode *prev = curr->prev;
        while (NULL != prev)
        {
          if ( (curr->time <= time) && (time <= prev->time) )
          {
            //Serial.println("+");
            return;
          }
          curr = prev;
          prev = curr->prev;
        }
        // shouldn't happen, the first node is too old!!
        memcpy(node->fields, curr->fields, sizeof(memoryNode::fields));
        node->quality = curr->quality;
        //Serial.println("<");
      }
    }

  private:
    memoryList(const memoryList &);
    memoryList&operator=(const memoryList &);

    memoryPool  cache;
    memoryNode *nodes;

};

static memoryList memory;


static void Memory_setup()
{
  // compute #nodes
  unsigned long num_nodes = 1 + (unsigned long)( ceilf(PROBE_DELAY * PROBE_FREQUENCY) );
  if (num_nodes <= 2)
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
  memory_last_time = Medium_GetCurrentTime();
}

// memory priming function
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

static void Memory_loop()
{
  const float t = Medium_GetCurrentTime();
  if (t - memory_last_time >= memory_delta_t)
  {
    memory.update();
    memory_last_time = t;
  }
}

static void Servo_loop()
{
  const float t = Medium_GetCurrentTime();
  memoryNode  node;
  memory.interpolate(t, &node);
  servo.write( 180.0f * Medium::Triangle(t, 10.0f) );
  servo_last_time = t;
}

// global setup
void setup()
{
  Serial.begin(baudRate);
  while (!Serial);
  Serial.print(F("sizeof(memoryNode) = ")); Serial.println(sizeof(memoryNode));
  Serial.print(F("sizeof(memoryNode::fields) = ")); Serial.println(sizeof(memoryNode::fields));

  Memory_setup();
  Servo_attach();
}

// input processing
void processInput()
{
  const unsigned numWords = medium.splitInput();
  Serial.print(F("#words=")); Serial.println(numWords);
  medium.resetInput();
}

// global loop
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

// global I/O
void serialEvent()
{
  medium.serialEventCallback();
}


