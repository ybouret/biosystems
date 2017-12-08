//_____________________________________________________________________________
//
// include required headers
//_____________________________________________________________________________
#include <Servo.h>

//_____________________________________________________________________________
//
//
// Board Setup, and servo with its init angle
//
//_____________________________________________________________________________
static const int            pinServo        = 9;      //!< Servo command channel (output PWM)
static const int            pinFrame        = 12;     //!< Frame capture control (output PWM)
static const int            pinFn           = 1;      //!< Normal   Force    (input Analog)
static const int            pinFt           = 4;      //!< Tangentil Force   (input Analog)
static const unsigned long  baudrate        = 115200; //!< check with Serial monitor/Python

static Servo                servo;
static const float          servo_angle_init  = 90.0f; //!< resting angle

//_____________________________________________________________________________
//
//
// Global Functions and Macros
//
//_____________________________________________________________________________

//-----------------------------------------------------------------------------
// DELAY: the delay time in seconds
// NODES: the number of nodes to reach the delay => sampling precision
//        is about DELAY/NODES
//-----------------------------------------------------------------------------
#define DELAY      1
#define NODES      100

//-----------------------------------------------------------------------------
//
// timing macros, GetCurrentTime is in seconds
//
//-----------------------------------------------------------------------------
#define TSYS()          (micros())
#define TSYS2TIME(tmx)  ( 1.0e-6f * (float)(tmx))
#define GetCurrentTime() TSYS2TIME( TSYS( ) )

//-----------------------------------------------------------------------------
//
// Forces: supposed to be between 0 and 5 volts
//
//-----------------------------------------------------------------------------
static float FRange  = 5.0f;  //!< force range in Volts
static float FnInit  = 2.50f; //!< initial normal force reading (means 0 force)
static float FtInit  = 2.50f; //!< initial tangential force reading (means 0 force)



//_____________________________________________________________________________
//
//
// Storage Management
//
//_____________________________________________________________________________

//-----------------------------------------------------------------------------
// we compute the store_rate so that
// (NODES-1)*store_rate > store_delay, resolution is microseconds
//-----------------------------------------------------------------------------
static const float         store_delay                   = (float)(DELAY);
static const unsigned long store_rate_extra_microseconds = 100;
static const unsigned long store_rate_microseconds       = (unsigned long)( (1.0e6f * store_delay) / (NODES - 1) + store_rate_extra_microseconds);
static const float         store_rate                    = ( (float)store_rate_microseconds) * 1.0e-6f;
static float               store_last_time               = 0.0f;

// doubly linked node
struct Node
{
  struct Node *next;  //!< for list
  struct Node *prev;  //!< for list
  float        time;  //!< last acquisition time
  float        angle; //!< last read angle
  float        Fn;    //!< normal force
  float        Ft;    //!< tangential force
};

// doubly linked list
struct List
{
  struct Node  *head;
  struct Node  *tail;
  unsigned      size;
};

//! device memory dependent #NODES
static struct Node nodes[NODES];

//! the data store
static struct List store = { NULL, NULL, 0 };

//-----------------------------------------------------------------------------
// used to initialize the store
//-----------------------------------------------------------------------------
static inline void __store_init(struct Node *node,
                                const float  time)
{
  node->prev  = node->next = 0;
  node->time  = time;
  node->angle = servo_angle_init;
  node->Fn    = FnInit;
  node->Ft    = FtInit;
  if (store.size <= 0)
  {
    store.head = store.tail = node;
  }
  else
  {
    store.tail->next = node;
    node->prev       = store.tail;
    store.tail       = node;
  }
  ++store.size;
}

//-----------------------------------------------------------------------------
// replace the last value, and put it in front of the list
//-----------------------------------------------------------------------------
static inline void StoreUpdate(const float time,
                               const float angle,
                               const float Fn,
                               const float Ft)
{
  struct Node *node = store.tail;

  // unlink tail
  store.tail       = node->prev;
  store.tail->next = NULL;
  node->prev       = NULL;

  // push_front the node
  node->next       = store.head;
  store.head->prev = node;
  store.head       = node;

  // update node data
  node->time   = time;
  node->angle  = angle;
  node->Fn     = Fn;
  node->Ft     = Ft;
}

//-----------------------------------------------------------------------------
// used to debug only
//-----------------------------------------------------------------------------
static inline void list_print()
{
  Serial.print("#");    Serial.print(store.size);
  const float t = GetCurrentTime();
  Serial.print(" t="); Serial.print(t);
  Serial.print("->"); Serial.print(t - store_delay);

  for (const struct Node *node = store.head; node != NULL; node = node->next)
  {
    Serial.print(" (");
    Serial.print(node->time);
    //Serial.print(",");Serial.print(node->angle);
    Serial.print(","); Serial.print(node->Fn);
    Serial.print(")");
  }
  Serial.println("");
}

//-----------------------------------------------------------------------------
//
//! initialize the store with "resting" nodes
//
//-----------------------------------------------------------------------------
static inline void StoreSetup()
{
  for (int i = 0; i < NODES; ++i)
  {
    const float time  = -((i + 1) * store_rate);
    __store_init( &nodes[i], time );
  }

}

//-----------------------------------------------------------------------------
//
//! read the forcing in Volts: range is 0:FnRange
//
//-----------------------------------------------------------------------------
static const float AnalogInputToVolts = (FRange / 1023.0f);
#define READ_FORCE_ON(PIN) ( AnalogInputToVolts * ( (float) analogRead(PIN) ) )
#define READ_FORCE_N() READ_FORCE_ON(pinFn)
#define READ_FORCE_T() READ_FORCE_ON(pinFt)

//-----------------------------------------------------------------------------
// called during the loop() function: store time, angle, values
//-----------------------------------------------------------------------------
static inline void StoreLoop()
{
  static const float InputToVolts = 5.0f / 1023.0f;

  const float local_time = GetCurrentTime();
  if (local_time - store_last_time >= store_rate)
  {
    const float Fn = READ_FORCE_N();
    const float Ft = READ_FORCE_T();
    StoreUpdate(local_time, servo.read(), Fn, Ft);
    store_last_time = GetCurrentTime();
  }
}

//-----------------------------------------------------------------------------
// compute an interpolated node:time/angle/value
//-----------------------------------------------------------------------------
static inline void StoreQuery(struct Node *data)
{
  const float  t          = GetCurrentTime() - store_delay;
  const struct Node *curr = store.tail;

  data->time = t; //!< the local time is saved
  if (curr->time >= t)
  {
    // timing is too short, shouldn't happen: this is a failsafe
    data->angle = curr->angle;
    data->Fn    = curr->Fn;
    data->Ft    = curr->Ft;
    return;
  }
  else
  {
    // we bracket the value
    const struct Node *prev = curr->prev;
    while (1)
    {
      if (NULL == prev)
      {
        // couldn't bracket, delay is probably too small
        // this is a failsafe
        data->angle = curr->angle;
        data->Fn    = curr->Fn;
        data->Ft    = curr->Ft;
        return;
      }
      else
      {
        // test bracketing
        const float p_time = prev->time;
        const float c_time = curr->time;
        if (p_time >= t && t >= c_time)
        {
          // bracketed!
          const float dt_num  = t - p_time;
          const float dt_den  = c_time - p_time;
          data->Fn    = prev->Fn    + (dt_num * (curr->Fn    - prev->Fn   )) / dt_den;
          data->Ft    = prev->Ft    + (dt_num * (curr->Ft    - prev->Ft   )) / dt_den;
          data->angle = prev->angle + (dt_num * (curr->angle - prev->angle)) / dt_den;
          return;
        }
        else
        {
          curr = prev;
          prev = prev->prev;
        }
      }
    }
  }
}

//-----------------------------------------------------------------------------
// example of function
//-----------------------------------------------------------------------------
static const float the_coeff = 50.0f;
static inline float ThetaDot()
{
  struct Node data;
  StoreQuery(&data);
  return -the_coeff * (data.Fn - FnInit);
}

//-----------------------------------------------------------------------------
// suivi non retarde
//-----------------------------------------------------------------------------
const float alpha = 5.0;

static inline float ThetaSuivi()
{
  struct Node data;
  StoreQuery(&data);
  const float force_factor = (data.Fn - FnInit) / FRange; //-0.5:0.5
  const float shift_factor = 180.0f * force_factor;   //-90:90
  return  alpha * shift_factor + 90; //avec gbf
}

//_____________________________________________________________________________
//
//
// Servo command
//
//_____________________________________________________________________________
float servo_last_time   = 0.0f;
float servo_rate        = 0.00f; //!< servo interaction rate

//-----------------------------------------------------------------------------
// shortcut...
//-----------------------------------------------------------------------------
static inline void ServoRest()
{
  servo.write(servo_angle_init);
}

//-----------------------------------------------------------------------------
// Servo and its data initialisation
//-----------------------------------------------------------------------------
static inline void ServoSetup()
{
  ServoRest();
  servo_last_time = GetCurrentTime();
}

//-----------------------------------------------------------------------------
//
// Servo loop function, compute the new angle
//
//-----------------------------------------------------------------------------
static float old_t = 0.0f;
static float dt    = 0.0f;
static float theta = servo_angle_init;
static const float theta_min =   0.0f;
static const float theta_max = 180.0f;

static inline void ServoLoop()
{
  //__________________________________________________
  //
  // get timings
  //__________________________________________________
  const float local_time = GetCurrentTime();
  dt    = local_time - old_t;
  old_t = local_time;

  //__________________________________________________
  //
  // where theta is computed, AT EVERY CALL
  //__________________________________________________
  theta = ThetaSuivi();

  //__________________________________________________
  //
  // clamp theta value
  //__________________________________________________
  if (theta >= theta_max) theta = theta_max;
  if (theta <= theta_min) theta = theta_min;

  //__________________________________________________
  //
  // when servo is set, at every servo_rate
  //__________________________________________________

  if ( local_time - servo_last_time >= servo_rate )
  {
    servo.write( theta  );

    servo_last_time  = local_time;
  }

}

//_____________________________________________________________________________
//
//
// Serial IO
//
//_____________________________________________________________________________
static float serial_last_time = 0.0f;
static float serial_rate      = 0.1f; //!< in seconds

void SerialSetup()
{
  serial_last_time = GetCurrentTime();
}

void SerialLoop()
{
  const float local_time = GetCurrentTime();
  if ( local_time - serial_last_time >= serial_rate )
  {
    // print useful info on serial: always the time
    Serial.print("t="); Serial.print(local_time);

    // physics
    Serial.print(" alpha="); Serial.print( alpha );
    Serial.print(" theta="); Serial.print( servo.read()   );
    Serial.print(" Fn=");    Serial.print( READ_FORCE_N() );
    Serial.print(" Ft=");    Serial.print( READ_FORCE_T() );

    // done !
    Serial.println("");
    // update
    serial_last_time = local_time;
  }
}


//_____________________________________________________________________________
//
//
// Setup
//
//_____________________________________________________________________________
void setup()
{
  //_________________________________________________________________________
  //
  // Serial communication setup
  //_________________________________________________________________________
  Serial.begin(baudrate);

  //_________________________________________________________________________
  //
  // Servo communication setup: TODO check init values 900,2100
  //_________________________________________________________________________
  servo.attach(pinServo, 900, 2100);
  //servo.attach(pinServo, 600, 2250);

  //_________________________________________________________________________
  //
  // set Frame Capture to LOW
  //_________________________________________________________________________
  pinMode(pinFrame, OUTPUT);
  digitalWrite(pinFrame, LOW);

  //_________________________________________________________________________
  //
  // prepare the objects
  //_________________________________________________________________________
  StoreSetup();
  ServoSetup();
  SerialSetup();
  old_t = GetCurrentTime();
}


//_____________________________________________________________________________
//
// Main loop
//
//_____________________________________________________________________________
void loop()
{
#if 1
  StoreLoop();
  ServoLoop();
  SerialLoop();
#else
  const float sweep_deg = 500.0f * GetCurrentTime();
  const float sweep_rad = ((float)M_PI) * sweep_deg / 180.0f;
  const float angle     = 90.0f + 90.0f * sin( sweep_rad );
  Serial.println(angle);
  servo.write(angle);
  delay(10);
  //Serial.println(servo.read());
  SerialLoop();
#endif

}

