#include <Servo.h>
#include <medium.h>

static const int   pinServo   = 9;
static Servo       servo;
static const unsigned long   baudRate   = 115200;

long lo = 750;
long hi = 2250;
float period = 12.0f;
int   motion = 0;

Medium medium;

void setup()
{
  Serial.begin(baudRate);
  servo.attach(pinServo, lo, hi);
  servo.write(90.0f);
}

void processInput()
{
  const unsigned numWords = medium.splitInput();
  Serial.print("#words="); Serial.println(numWords);
  for (unsigned i = 0; i < numWords; ++i)
  {
    Serial.print("\t"); Serial.println(medium[i]);
  }

  if (2 == numWords)
  {
    const char *msg = medium[0];
    const char *arg = medium[1];

    if ( 0 == strcmp(msg, "period") )
    {
      const float value = atof(arg);
      if (value > 0)
      {
        period = value;
      }
      goto END;
    }

    if ( 0 == strcmp(msg, "motion" ) )
    {
      motion = atoi(arg);
      goto END;
    }
  }

END:
  // cleanup for next input
  medium.resetInput();
}

void loop()
{
  float t = Medium_GetCurrentTime();
  switch (motion)
  {
    case 0:
      servo.write( 180.0f * Medium::Triangle( t, period ) );
      break;

    case 1:
      servo.write( 90.0f + 90.0f * Medium::SineWave( t, period ) );
      break;

    case 2: {
        servo.write(0);
        delay(1000);
        t = Medium_GetCurrentTime();
        servo.write(180);
        while (Medium_GetCurrentTime() - t < 0.5)
        {
          Serial.print(servo.read());
          Serial.print(" ");
        }
        Serial.println(" ");
        motion = 0;
      } break;
  }

  if (medium.inputCompleted() )
  {
    processInput();
  }
}

void serialEvent()
{
  medium.serialEventCallback();
}

