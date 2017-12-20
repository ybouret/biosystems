#include <Servo.h>
#include <medium.h>

// Board setup
static const int servo_pin = 9; //!< Servo            | PWM input channel
static const int pinFn     = 1; //!< Normal Force     | analog input
static const int pinFt     = 4; //!< Tangential Force | analog input

static const unsigned long baudRate = 115200; //!< COM rate

// Medium setup
static Medium medium;

// Servo setup
static Servo servo;
static float servo_angle_mini = 0.0f;
static float servo_angle_maxi = 180.0f;
static float servo_angle_init = 90.0f;
static unsigned long servo_pwm_lower = 750;
static unsigned long servo_pwm_upper = 2250;

static void Servo_attach()
{
  servo.detach();
  servo.attach(servo_pin, servo_pwm_lower, servo_pwm_upper);
}

// global setup
void setup()
{
  Serial.begin(baudRate);
  Servo_attach();
}

// input processing
void processInput()
{
  const unsigned numWords = medium.splitInput();
  Serial.print("#words=");Serial.println(numWords);
  medium.resetInput();
}

// global loop
void loop()
{
  const float t = Medium_GetCurrentTime();
  servo.write( 180.0f * Medium::Triangle(t, 10));

  // I/O
  if (medium.inputCompleted() )
  {
    processInput();
  }
}

// global I/O
void serialEvent()
{
  medium.serialEventCallback();
}


