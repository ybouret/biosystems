#include <Servo.h>

static const int   pinServo   = 9;
static Servo       servo;
static const float angle_init = 90.0f;
static const unsigned long   baudrate   = 115200;

static float period     = 6.0f;

//-----------------------------------------------------------------------------
//
// timing macros, GetCurrentTime is in seconds
//
//-----------------------------------------------------------------------------
#define TSYS()          (micros())
#define TSYS2TIME(tmx)  ( 1.0e-6f * (float)(tmx))
#define GetCurrentTime() TSYS2TIME( TSYS( ) )

String  inputString = "";         // a String to hold incoming data
boolean stringComplete = false;  // whether the string is complete

#if 0
long lo = 544;
long hi = 2400;
#endif

long lo = 750;
long hi = 2250;


float amin = 0.0;
float amax = 180.0;

void setup()
{
  Serial.begin(baudrate);
  servo.attach(pinServo, lo, hi);
  servo.write(angle_init);
  inputString.reserve(200);
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

//float the_angle = 90.0f;

void process_message()
{
  static const char sep[] = " \t\r\n";
  char *data = (char *)inputString.c_str();
  char *msg  = strtok(data, sep);
  if (!msg)
    goto SHOW;
  else
  {
    char *arg  = strtok(NULL, sep);
    if (!arg) return;
    Serial.print("msg='"); Serial.print(msg); Serial.print("', arg='"); Serial.print(arg); Serial.println("'");

    if ( 0 == strcmp(msg, "period") )
    {
      const float tmp = atof(arg);
      if (tmp > 0)
      {
        period = tmp;
      }
      return;
    }

    if ( 0 == strcmp(msg, "amin") )
    {
      const float tmp = atof(arg);
      if (tmp >= 0.0f && tmp <= amax)
      {
        amin = tmp;
      }
      goto SHOW;
    }

    if ( 0 == strcmp(msg, "amax") )
    {
      const float tmp = atof(arg);
      if (tmp <= 180.0f && tmp >= amin)
      {
        amax = tmp;
      }
      goto SHOW;
    }

    if ( 0 == strcmp(msg, "lo") )
    {
      const  long tmp = atol(arg);
      if (tmp >= 0 && tmp <= hi)
      {
        lo = tmp;
        servo.detach();
        servo.attach(pinServo, lo, hi);
      }
      goto SHOW;
    }

    if ( 0 == strcmp(msg, "hi") )
    {
      const  long tmp = atol(arg);
      if (tmp >= 0 && tmp >= lo)
      {
        hi = tmp;
        servo.detach();
        servo.attach(pinServo, lo, hi);
      }
      goto SHOW;
    }
  }
SHOW:
  Serial.print("lo     = "); Serial.println(lo);
  Serial.print("hi     = "); Serial.println(hi);
  Serial.print("period = "); Serial.println(period);
  Serial.print("amin   = "); Serial.println(amin);
  Serial.print("amax   = "); Serial.println(amax);
}

float last_t = 0;
void loop()
{
  const float t = GetCurrentTime();

  {
    const float a = amin + (amax-amin)* Triangle(t);
    if ( t - last_t >= 0.5 )
    {
      last_t = t;
    }
    servo.write(a);
  }
  
  // print the string when a newline arrives:
  if (stringComplete) {
    Serial.println(inputString);
    process_message();
    // clear the string:
    inputString    = "";
    stringComplete = false;
  }

  //servo.write(the_angle);
}

/*
  SerialEvent occurs whenever a new data comes in the hardware serial RX. This
  routine is run between each time loop() runs, so using delay inside loop can
  delay response. Multiple bytes of data may be available.
*/
void serialEvent() {
  while (Serial.available()) {
    // get the new byte:
    char inChar = (char)Serial.read();
    // add it to the inputString:
    inputString += inChar;
    // if the incoming character is a newline, set a flag so the main loop can
    // do something about it:
    if (inChar == '\n') {
      stringComplete = true;
    }
  }
}


