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

void processInput()
{
  const unsigned numWords = medium.findInputWords();
  Serial.print("#words="); Serial.println(numWords);
  for(unsigned i=0;i<numWords;++i)
  {
      Serial.print("\t#"); Serial.print(i); Serial.print("=>'"); Serial.print( medium.getInputWord(i) ); Serial.println("'");
  }
  // cleanup for next input
  medium.resetInput();
}

void loop()
{
  const float t = Medium_GetCurrentTime();
  //servo.write( 180.0f * Medium::Triangle( t, period ) );
  servo.write( 90.0f + 90.0f * Medium::SineWave( t, period ) );
  if (medium.inputCompleted() )
  {
    Serial.print("received '"); Serial.print(medium.input()); Serial.println("'");
    processInput();
  }

}

void serialEvent()
{
  medium.processSerialInput();
}

