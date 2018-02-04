////////////////////////////////////////////////////////////////////////////////
//
// program to impose motion to fish and measure forces at the same time
//
////////////////////////////////////////////////////////////////////////////////

#include <Servo.h>
#include <medium.h>


//! OUTPUT_FREQUENCY in Hz
#define OUTPUT_FREQUENCY  2.0f

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
        Serial.print(F(" q="));      Serial.println(quality);
#endif
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



static void on_period(const char *value)
{
  const float tmp = atof(value);
  if (tmp > 0)
  {
    fish.period = tmp;
  }
  Serial.print(F("period=")); Serial.println(fish.period);
}

static void on_amplitude(const char *value)
{
  const float tmp = atof(value);
  if (tmp > 0)
  {
    fish.amplitude = tmp;
  }
  Serial.print(F("amplitude=")); Serial.println(fish.amplitude);
}

static void on_motion(const char *value)
{
  if (Medium_streq(value, "tri"))
  {
    fish.motion = Medium::TriangleWave; Serial.println("TriangleWave");
    return;
  }

  if (Medium_streq(value, "cos"))
  {
    fish.motion = Medium::CosWave;  Serial.println("CosWave");
    return;
  }

}

static const char PROGMEM periodID[] = "period";
static const Medium::Callback  PROGMEM periodCB   = on_period;

static const Medium::Parameter
//PROGMEM
__params[]  =
{
  MEDIUM_PARAM(period),
  MEDIUM_PARAM(motion),
  MEDIUM_PARAM(amplitude)
};

static const char PROGMEM string_a[] = "period";
static const char PROGMEM string_b[] = "motion";

struct Info
{
  const char * PROGMEM name;
  int                  indx;  
};

static const Info PROGMEM infos[] =
{
  { string_a, 67}
};

static const char * const PROGMEM Data [] =
{
  string_a,
  string_b
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
  servo.write(90.0f);

  // prepare the timings
  fish.last_output = Medium_GetCurrentTime();

#if 0
  Serial.print("sizeof(param)="); Serial.println(sizeof(Medium::Parameter));
  for (size_t i = 0; i < sizeof(__params) / sizeof(__params[0]); ++i)
  {
    Medium::Parameter p = {NULL, ""};
    //memcpy_P(&p, pgm_read_word( (&__params[i]) ), sizeof(Medium::Parameter));
    memcpy(&p, &__params[i], sizeof(Medium::Parameter));
    Serial.print("name="); Serial.println(p.name);
  }
#endif

#if 0
  for (size_t i = 0; i < sizeof(Data) / sizeof(Data[0]); ++i)
  {
    char buffer[16];
    strcpy_P(buffer, (const char *)pgm_read_word( Data + i ));
    Serial.print("Data="); Serial.println(buffer);
  }
#endif

 #if 0
  for (size_t i = 0; i < sizeof(infos) / sizeof(infos[0]); ++i)
  {
    char buffer[16];
    strcpy_P(buffer, (const char *)pgm_read_word( infos[i].name ));
    Serial.print("Data="); Serial.println(buffer);
  }
#endif
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
#if 0
  if ( medium.inputCompleted() )
  {
    Serial.println("---> got input");

    if (2 == medium.splitInput(NULL))
    {
      const char   *cmd  = medium.getField(0);
      const unsigned num_params = sizeof(__params) / sizeof(__params[0]);
      for (unsigned i = 0; i < num_params; ++i)
      {

        const Medium::Parameter &info = __params[i];
        Serial.print("\ttesting"); Serial.println(info.name);
        if ( Medium_streq(info.name, cmd) )
        {
          const char *value_string = medium.getField(1);
          info.proc(value_string);
          goto END_INPUT;
        }
      }
    }

END_INPUT:
    medium.resetInput();
  }
#endif
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



