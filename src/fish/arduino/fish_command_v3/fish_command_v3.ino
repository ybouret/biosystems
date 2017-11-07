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

float global_last_time = 0.0f; //!< last time in main loop
float global_rate      = 0.2f; //!< main loop processing rate

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
    // and rest a little...
    //_________________________________________________________________________
    delay(1000);

    global_last_time = GetCurrentTime();
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

    // check if we process something
    if( (local_curr_time-global_last_time)>=global_rate && local_curr_time <= 10 )
    {   
        Serial.print(local_curr_time) ;
        Serial.print("\t");
        Serial.print(global_last_time);
        Serial.println("");

        global_last_time = GetCurrentTime();
    }

}