#include <Servo.h>
#include <medium.h>

static const int   pinServo   = 9;
static Servo       servo;
static const unsigned long   baudRate   = 115200;

long lo = 750;
long hi = 2250;
float period = 10.0f;

Medium medium;

void setup()
{
  Serial.begin(baudRate);
  servo.attach(pinServo, lo, hi);
  servo.write(90.0f);
}

void p(const char *fmt, ... ) {
  char buf[128]; // resulting string limited to 128 chars
  va_list args;
  va_start (args, fmt );
  vsnprintf(buf, 128, fmt, args);
  va_end (args);
  Serial.print(buf);
}

void processInput()
{
  const unsigned numWords = medium.findInputWords();
  medium.print("#words = %4u\n", numWords);
  /*
  p("#words = %4u\n", numWords);
  for (unsigned i = 0; i < numWords; ++i)
  {
    Serial.print("\t#"); Serial.print(i); Serial.print("=>'"); Serial.print( medium.getInputWord(i) ); Serial.println("'");
  }
  */
  // cleanup for next input
  medium.resetInput();
}

void loop()
{
  const float t = Medium_GetCurrentTime();
  servo.write( 180.0f * Medium::Triangle( t, period ) );
  //servo.write( 90.0f + 90.0f * Medium::SineWave( t, period ) );
  if (medium.inputCompleted() )
  {
    processInput();
  }

}

void serialEvent()
{
  medium.processSerialInput();
  medium.processSerialOutput();
}

