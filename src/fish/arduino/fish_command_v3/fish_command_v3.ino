//_____________________________________________________________________________
//
//
//_____________________________________________________________________________
#include <Servo.h>

//_____________________________________________________________________________
//
//
// Board Setup
//
//_____________________________________________________________________________
const int            pinServo        = 9;  //!< Servo command
const int            pinFrame        = 12; //!< Frame capture control
const int            pinVoltageInput = 2;  //!< input Voltage in 2
const unsigned long  baudrate        = 115200;

//_____________________________________________________________________________
//
//
// global functions and macros
//
//_____________________________________________________________________________
#define TSYS()          (micros())
#define TSYS2TIME(tmx)  ( 1.0e-6f * (float)(tmx))
#define GetCurrentTime() TSYS2TIME( TSYS( ) )

//_____________________________________________________________________________
//
//
// Servo/fish commamnd
//
//_____________________________________________________________________________
Servo servo;
float servo_angle_init  = 90.0f; //!< resting angle
float servo_angle_shift =  0.0f; //!< adjust in real world
float servo_swim_time   =  0.0f; //!< time where it starts to swim
float servo_last_time   =  0.0f;
float servo_last_angle  = servo_angle_init;

float servo_rate       = 0.2f; //!< servo interaction rate


static inline void ServoSetAngle(const float angle)
{
    servo.write(angle+servo_angle_shift);
}

static inline void ServoRest()
{
    ServoSetAngle(servo_angle_init);
}

static inline void ServoSetup()
{
  ServoRest();
  servo_last_time = GetCurrentTime();  
}

static inline void ServoLoop()
{
    const float local_time = GetCurrentTime();
    if( local_time - servo_last_time >= servo_rate )
    {
      Serial.print("servo@t=");Serial.println(local_time);
      servo_last_time = GetCurrentTime();
    }
}



//_____________________________________________________________________________
//
//
// chained list management
//
//_____________________________________________________________________________
#define NODES             8
static const float store_delay = 1.0f;
static const float store_rate  = store_delay/NODES;
static float       store_last_time   = 0.0f; //!< 

struct Node 
{
  struct Node *next;
  struct Node *prev;
  float        time;
  float        angle;
  float        value;
};

struct List 
{
  struct Node  *head;
  struct Node  *tail;
  unsigned      size;
};

struct Node nodes[NODES];
struct List history = { NULL, NULL, 0 };

static inline void list_push_back(struct Node *node, const float time)
{
   node->prev  = node->next = 0;
   node->time  = time;
   node->angle = servo_angle_init;
   node->value = 0;
   if(history.size<=0)
   {
      history.head = history.tail = node;
   }
   else 
   {
      history.tail->next = node;
      node->prev         = history.tail;
      history.tail       = node;
   }
   ++history.size;
}

static inline void list_store(const float time, const float angle, const float value)
{
    struct Node *node = history.tail;
    
    // unlink tail
    history.tail       = node->prev;
    history.tail->next = NULL;
    node->prev = NULL;
   
    // push front the node
    node->next         = history.head;
    history.head->prev = node;
    history.head = node;
    node->time  = time;
    node->angle = angle;
    node->value = value;
}

static inline void list_print()
{
  Serial.print("#"); Serial.print(history.size);
  for(const struct Node *node = history.head;node!=NULL;node=node->next)
  {
    Serial.print(" (");
    Serial.print(node->time);
    Serial.print(",");Serial.print(node->angle);
    Serial.print(",");Serial.print(node->value);
    Serial.print(")");
  }
  Serial.println("");
}

static inline void StoreSetup()
{
    for(int i=0;i<NODES;++i)
    {
        list_push_back( &nodes[i], -(i*store_rate)  );
        history.tail->angle = servo_angle_init + 45 * sin( 2*3.14 * history.tail->time / store_delay );
    }
    
}

static inline void StoreLoop()
{
    const float local_time = GetCurrentTime();
    if(local_time-store_last_time>=store_rate)
    {
      list_store(local_time,servo.read(),analogRead(pinVoltageInput));
      list_print();
      store_last_time = GetCurrentTime();
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
    servo.attach(pinServo,900,2100);
    

    //_________________________________________________________________________
    //
    // set Frame Capture to LOW
    //_________________________________________________________________________
    pinMode(pinFrame,OUTPUT);
    digitalWrite(pinFrame,LOW);

    //_________________________________________________________________________
    //
    // prepare the objects
    //_________________________________________________________________________
    StoreSetup();
    ServoSetup();
}


//_____________________________________________________________________________
//
//
//
//_____________________________________________________________________________
void loop()
{
    StoreLoop();
    ServoLoop();
}

