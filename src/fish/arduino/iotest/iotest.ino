#include <Servo.h>
#include <medium.h>

static const int   pinServo   = 9;
static Servo       servo;
static const unsigned long   baudRate   = 115200;

long lo = 750;
long hi = 2250;
float period = 12.0f;

Medium medium;

void setup()
{
  Serial.begin(baudRate);
  servo.attach(pinServo, lo, hi);
  servo.write(90.0f);
  Serial.print("sizeof(List)="); Serial.println(sizeof(Medium::List< Medium::NodeOf<int> >));
  Serial.print("sizeof(CharNode)="); Serial.println(sizeof(Medium::NodeOf<char>));
}

void processInput()
{
  const unsigned numWords = medium.findInputWords();
  medium.print("#words = %4u\n", numWords);
  for(unsigned i=0;i<numWords;++i)
  {
    medium.print("\t=> %s\n", medium.getInputWord(i) );
  }

  if(2==numWords)
  {
    const char *msg = medium.getInputWord(0);
    const char *arg = medium.getInputWord(1);
    
  }

  END:
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
    processInput();
  }
  medium.loopCallback();

}

void serialEvent()
{
  medium.serialEventCallback();
}

