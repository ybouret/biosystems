////////////////////////////////////////////////////////////////////////////////
//
// program to impose motion to fish and measure forces at the same time
//
////////////////////////////////////////////////////////////////////////////////

#include <Servo.h>
#include <medium.h>


//! OUTPUT_FREQUENCY in Hz
#define OUTPUT_FREQUENCY  5.0f

////////////////////////////////////////////////////////////////////////////////
//
// Board setup
//
////////////////////////////////////////////////////////////////////////////////

#define pinServo  9
#define pinFn     1
#define pinFt     4
#define baudRate  115200

////////////////////////////////////////////////////////////////////////////////
//
// I/O data
//
////////////////////////////////////////////////////////////////////////////////


// FORCES_VOLTAGE in Volts
#define FORCES_VOLTAGE 5.0f

class RoboFish : public Servo
{
  public:
    typedef float (*SwimMotion)(float, const float);
    const Medium::AnalogReader read_force;

    const float dt_output;   //!< 1/OUTPUT_FREQUENCY
    float       last_output; //!< time to write data
    float       period;      //!< motion frequency
    float       amplitude;   //!< motion amplitude
    SwimMotion  motion;      //!< default: Medium::Triangle

    //! constructor: initialize data
    inline explicit RoboFish() throw() :
      read_force( FORCES_VOLTAGE ),
      dt_output(1.0f / OUTPUT_FREQUENCY),
      last_output(-1),
      period(5.0f),
      amplitude(90.0f),
      motion( Medium::Triangle )
    {}

    //! destructor
    inline virtual ~RoboFish() throw() {
      /* nothing to do*/
    }


    inline void loop()
    {
      const float t = Medium_GetCurrentTime();

      // set the servo
      const float angle = 90.0f + amplitude * motion(t, period);
      write(angle);

      // probe the forces
      const float real_dt = t - last_output;
      if (real_dt >= dt_output)
      {
        const float Fn = read_force(pinFn);
        const float Ft = read_force(pinFt);
        Serial.print(F("t="));      Serial.print(t);
        Serial.print(F(" angle=")); Serial.print(angle);
        Serial.print(F(" Fn="));    Serial.print(Fn);
        Serial.print(F(" Ft="));    Serial.print(Ft);
        Serial.println(F(""));
        last_output = t;
      }
    }

  private:
    MEDIUM_DISABLE_COPY_AND_ASSIGN(RoboFish);
};

static RoboFish  fish;
static Medium    medium;

////////////////////////////////////////////////////////////////////////////////
//
// input processing
//
////////////////////////////////////////////////////////////////////////////////

//! table of commands into flash memory
static  const char * const PROGMEM commands[]  =
{
  "period"
};

//! command aliases
#define __PERIOD commands[0]

static void processInput()
{
  const unsigned numWords = medium.splitInput();
  if (2 == numWords)
  {
    const char *cmd = medium[0];
    if ( Medium_streq(__PERIOD, cmd) )
    {
      const float tmp = atof(medium[1]);
      if (tmp > 0)
      {
        fish.period = tmp;
      }
      goto END_INPUT;
    }
  }

  goto END_INPUT;
END_INPUT:
  medium.resetInput();
}

////////////////////////////////////////////////////////////////////////////////
//
// GLOBAL setup
//
////////////////////////////////////////////////////////////////////////////////
void setup()
{
  // the proper way to initialize serial I/O initialize, see doc
  Serial.begin(baudRate);
  while (!Serial)
    ;

  // prepare the servo
  Servo &servo = fish;
  servo.detach();
  servo.attach(pinServo); // TODO: check PWN parameters
  servo.write(90.0f);

  // prepare the timings
  fish.last_output = Medium_GetCurrentTime();
}

////////////////////////////////////////////////////////////////////////////////
//
// GLOBAL loop
//
////////////////////////////////////////////////////////////////////////////////
void loop()
{
  fish.loop();
  // I/O
  if (medium.inputCompleted() )
  {
    processInput();
  }
}

////////////////////////////////////////////////////////////////////////////////
//
// global I/O
//
////////////////////////////////////////////////////////////////////////////////
void serialEvent()
{
  medium.serialEventCallback();
}



