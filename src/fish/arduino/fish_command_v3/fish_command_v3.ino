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
const int            pinServo      = 9;  //!< Servo command
const int            pinFrame      = 12; //!< Frame capture control
const unsigned long  baudrate      = 115200;

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
// chained list management
//
//_____________________________________________________________________________
float store_last_time  = 0.0f; //!< 
float store_rate       = 0.5f; //!< storing new position

struct Node 
{
  struct Node *next;
  struct Node *prev;
  float        time;
  float        angle;
};

struct List 
{
  struct Node  *head;
  struct Node  *tail;
  unsigned      size;
};

#define NODES 8
struct Node nodes[NODES];
struct List history = { NULL, NULL, 0 };

static void list_push_back(struct Node *node)
{
   node->prev  = node->next = 0;
   node->time  = 0;
   node->angle = servo_angle_init;
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

static void list_store(const float time, const float angle)
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
}

static void list_print()
{
  Serial.print("#"); Serial.print(history.size);
  for(const struct Node *node = history.head;node!=NULL;node=node->next)
  {
    Serial.print(" (");Serial.print(node->time);Serial.print(",");Serial.print(node->angle);Serial.print(")");
  }
  Serial.println("");
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
    // Serial communication setup: TODO check init values 900,2100
    //_________________________________________________________________________
    servo.attach(pinServo,900,2100);
    ServoSetAngle(servo_angle_init);

    //_________________________________________________________________________
    //
    // set Frame Capture to LOW
    //_________________________________________________________________________
    pinMode(pinFrame,OUTPUT);
    digitalWrite(pinFrame,LOW);

    //_________________________________________________________________________
    //
    // prepare the list
    //_________________________________________________________________________
    for(unsigned i=0;i<NODES;++i)
    {
        list_push_back( &nodes[i] );
    }

    //_________________________________________________________________________
    //
    // and rest a little...
    //_________________________________________________________________________
    delay(1000);

    servo_last_time = GetCurrentTime();
    store_last_time = GetCurrentTime();
}


//_____________________________________________________________________________
//
//
//
//_____________________________________________________________________________
void loop()
{
    // get current time
    const float local_curr_time = GetCurrentTime();

    // check servo
    if( (local_curr_time-servo_last_time)>=servo_rate && local_curr_time <= 10 )
    {   
        Serial.print("Current Time:");
        Serial.print(local_curr_time) ;
        Serial.println("");
        servo_last_time = GetCurrentTime();
    }

    // check store
    if( (local_curr_time-store_last_time)>=store_rate )
    {
        Serial.print("Store ");
        Serial.print(local_curr_time);
        Serial.print(",");
        Serial.print(servo_last_angle);
        Serial.println("");
        list_store(local_curr_time,servo_last_angle);
        store_last_time = local_curr_time;
        list_print();
    }
  
    

}
