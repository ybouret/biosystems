////////////////////////////////////////////////////////////////////////////////
//
// program to impose motion to fish and measure forces at the same time
//
////////////////////////////////////////////////////////////////////////////////

#include <Servo.h>
#include <medium.h>


//! OUTPUT_FREQUENCY in Hz
#define OUTPUT_FREQUENCY  10.0f

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
      motion( Medium::CosWave )
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
        const float Fn      = read_force(pinFn);
        const float Ft      = read_force(pinFt);
        const int   quality = (int)floorf( (100.0f * dt_output / real_dt) + 0.5f );

#if 1
        Serial.print(F("t="));      Serial.print(t);
        Serial.print(F(" angle=")); Serial.print(angle);
        Serial.print(F(" Fn="));    Serial.print(Fn);
        Serial.print(F(" Ft="));    Serial.print(Ft);
        Serial.print(F(" q="));     Serial.println(quality);
#endif
        last_output = t;
        //Serial.println(F("#comment"));
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


static void on_period(const char *value)
{
  const float tmp = atof(value);
  if (tmp > 0)
  {
    fish.period = tmp;
  }
  Serial.print(F("#period=")); Serial.println(fish.period);
}

static void on_amplitude(const char *value)
{
  const float tmp = atof(value);
  if (tmp > 0)
  {
    fish.amplitude = tmp;
  }
  Serial.print(F("#amplitude=")); Serial.println(fish.amplitude);
}

static const char PROGMEM __tri[] = "tri";
static const char PROGMEM __cos[] = "cos";

static void on_motion(const char *value)
{

  if (Medium_streq_P(value, __tri))
  {
    fish.motion = Medium::TriangleWave; Serial.println(F("#TriangleWave"));
    return;
  }

  if (Medium_streq_P(value, __cos))
  {
    fish.motion = Medium::CosWave;  Serial.println(F("#CosWave"));
    return;
  }

}

static const Medium::Parameter __params[]  =
{
  MEDIUM_PARAM(period),
  MEDIUM_PARAM(motion),
  MEDIUM_PARAM(amplitude)
};


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
  servo.write(90.0f+fish.amplitude);

  // prepare the timings
  fish.last_output = Medium_GetCurrentTime();

  //send a comment line to help Serial sync!!!
  Serial.println(F("#ready!"));
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
  medium.processInput(MEDIUM_PARAMETERS(__params));
}

////////////////////////////////////////////////////////////////////////////////
//
// Global I/O
//
////////////////////////////////////////////////////////////////////////////////
void serialEvent()
{
  medium.serialEventCallback();
}



