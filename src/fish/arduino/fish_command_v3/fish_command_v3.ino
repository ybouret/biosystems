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
static const int            pinValueInput   = 0;      //!< Value Input a value   (input Analog)
static const unsigned long  baudrate        = 115200; //!< check with Serial monitor/Python

static Servo                servo;
static const float          servo_angle_init  = 90.0f; //!< resting angle

//_____________________________________________________________________________
//
//
// global functions and macros
//
//_____________________________________________________________________________

//-----------------------------------------------------------------------------
// DELAY: the delay time..
// NODES: the number of nodes to reach the delay => sampling precision
//        is about DELAY/NODES
//-----------------------------------------------------------------------------
#define DELAY      1.5
#define NODES      50

//-----------------------------------------------------------------------------
// timing macros
//-----------------------------------------------------------------------------
#define TSYS()          (micros())
#define TSYS2TIME(tmx)  ( 1.0e-6f * (float)(tmx))
#define GetCurrentTime() TSYS2TIME( TSYS( ) )



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
  struct Node *next;
  struct Node *prev;
  float        time;
  float        angle;
  float        value;
};

// doubly linked list
struct List
{
  struct Node  *head;
  struct Node  *tail;
  unsigned      size;
};

// device memory dependent #NODES
static struct Node nodes[NODES];

// the data store
static struct List store = { NULL, NULL, 0 };

//-----------------------------------------------------------------------------
// used to initialize the store
//-----------------------------------------------------------------------------
static inline void store_push_back(struct Node *node, const float time, const float angle)
{
  node->prev  = node->next = 0;
  node->time  = time;
  node->angle = angle;
  node->value = 0;
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
static inline void list_store(const float time, const float angle, const float value)
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
  node->value  = value;
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
    Serial.print(","); Serial.print(node->value);
    Serial.print(")");
  }
  Serial.println("");
}

//-----------------------------------------------------------------------------
// initialize the store
//-----------------------------------------------------------------------------
static inline void StoreSetup()
{
  for (int i = 0; i < NODES; ++i)
  {
    const float time  = -((i + 1) * store_rate);
    store_push_back( &nodes[i], time, servo_angle_init  );
  }

}

//-----------------------------------------------------------------------------
// called during the loop() function: store time, angle, value
//-----------------------------------------------------------------------------
static inline void StoreLoop()
{
  const float local_time = GetCurrentTime();
  if (local_time - store_last_time >= store_rate)
  {
    const int    analogValue = analogRead(pinValueInput);        //!< in 0:1023
    const float  value       = ( (float)analogValue ) / 1023.0f; //!< in 0:1.0f
    list_store(local_time, servo.read(), value);
    //list_print();
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
    data->value = curr->value;
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
        // couldn't bracket, shouldn't happen: this is a failsafe
        data->angle = curr->angle;
        data->value = curr->value;
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
          const float p_value = prev->value;
          const float c_value = curr->value;
          const float dt_num  = t - p_time;
          const float dt_den  = c_time - p_time;
          data->value = prev->value + (dt_num * (curr->value - prev->value)) / dt_den;
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
static const float alpha = 50;
static inline float ThetaDot()
{
  struct Node data;
  StoreQuery(&data);
  return -alpha * (data.value - 0.5);
}


//_____________________________________________________________________________
//
//
// Servo command
//
//_____________________________________________________________________________
float servo_last_time   = 0.0f;
float servo_last_angle  = servo_angle_init;
float servo_rate        = 0.1f; //!< servo interaction rate

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
// Servo loop function, compute the new angle
//-----------------------------------------------------------------------------
static float old_t = 0.0;
static float theta = servo_angle_init;

static inline void ServoLoop()
{
  const float local_time = GetCurrentTime();

  {
    const float dt = local_time - old_t;
    old_t  = local_time;
    theta += ThetaDot() * dt;
    if (theta >= 180.0f) theta = 180.0f;
    if (theta <= 0.0f)   theta = 0.0f;
  }

  if ( local_time - servo_last_time >= servo_rate )
  {
    // do something with the servo
    Serial.println(theta);
    servo.write( theta  );

    // save info
    servo_last_time  = GetCurrentTime();
    servo_last_angle = theta;
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
  //servo.attach(pinServo,900,2100);
  servo.attach(pinServo, 600, 2250);

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
  old_t = GetCurrentTime();
}


//_____________________________________________________________________________
//
//
//
//_____________________________________________________________________________
void loop()
{
#if 1
  StoreLoop();
  ServoLoop();
#else
  servo.write(0);
  delay(500);
  servo.write(180);
  delay(500);
#endif

}

