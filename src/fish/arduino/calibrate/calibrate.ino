#include <Servo.h>

static const int   pinServo   = 9;
static Servo       servo;
static const float angle_init = 90.0f;
static const unsigned long   baudrate   = 57600;

static const float period     = 5.0f;

//-----------------------------------------------------------------------------
//
// timing macros, GetCurrentTime is in seconds
//
//-----------------------------------------------------------------------------
#define TSYS()          (micros())
#define TSYS2TIME(tmx)  ( 1.0e-6f * (float)(tmx))
#define GetCurrentTime() TSYS2TIME( TSYS( ) )

void setup()
{
  Serial.begin(baudrate);
  servo.attach(pinServo, 900, 2100);
  servo.write(angle_init);
}

// 0->1 on -T/2:0, 1->0 on 0:T/2
float Triangle(float t)
{
  const float half = 0.5f * period;
  t = t - period * floorf( t / period + 0.5f );
  if (t <= 0)
  {
    return (t + half) / half;
  }
  else
  {
    return (half - t) / half;
  }
}

float last_t = 0;
void loop()
{
  const float t = GetCurrentTime();
  if (t <= 20)
  {
    const float a = 180.0f * Triangle(t);
    if ( t - last_t >= 0.5 )
    {
      last_t = t;
      Serial.println("+");
    }
    servo.write(a);
  }
  else
  {
    servo.write(90);
  }

}

