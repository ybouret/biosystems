#include <Servo.h>
#include <math.h>

//______________________________________________________________________________
//
//
// communication variables
//
//______________________________________________________________________________
unsigned long     baudrate       = 115200;
#define           pin_servo        9

#if 0
const int         pin_servo      = 9;
int      pin_interFrame     = 12;
int      pin_forceSensor    = 1;
int      pin_tensionEntree  = 3;
int      pin_tensionSortie  = 2;
#endif



//______________________________________________________________________________
//
//
// Simulation parameters, user defined
//
//______________________________________________________________________________
Servo    myservo;
float    frequency          = 1.0;     // servo frequency, in Hz
float    amplitude          = 45.0;   // in degrees
unsigned num_periods        = 5;       // number of servo periods during motion
unsigned points_per_period  = 20;      // number of saved points per cycle
float    resting_time       = 0.5;     // resting initial time
unsigned points_during_rest = 10;      // points to save during resting time
float    angle_init         = 90.0f;
float    angle_shift        =  0.0f;

//______________________________________________________________________________
//
//
// Simulation parameters, automatic
//
//______________________________________________________________________________
float    dt_rest      = resting_time/points_during_rest;
unsigned num_iter     = points_during_rest + points_per_period * num_periods;
unsigned iter         = 0;

float period          = 1.0/frequency;
float dt_swim         = period/points_per_period;
const float omega     = 2.0*M_PI*frequency;
const float start_angle = angle_init+angle_shift;



unsigned long     tsys_start     = 0;      // in microseconds
unsigned long     tsys           = 0;
float             t              = 0;
float             t_last         = 0;
float             t_save         = 0.25;

bool              resting     = true;
float             t_swim      = 0;


#define TSYS()       (micros())
#define TSYS2TIME(T) (1e-6* ((float)(T)))


// write resting position
void myservo_rests()
{
    myservo.write(angle_init);
}

// write swimming position
void myservo_swims()
{
    // ellapsed time since swimming order
    const float t_run = t - t_swim;
    const float angle = amplitude * sin( omega * t_run ) + start_angle;
    myservo.write( (int) angle );
}

//______________________________________________________________________________
//
//
// Initialisation: probes, servo and Python info
//
//______________________________________________________________________________
void setup()
{
    // communication settings
    Serial.begin(baudrate);
    
    // initialize servo
    myservo.attach(pin_servo,900,2100);
    myservo_rests();
    delay(500);
    
    // initialize
    tsys_start = TSYS();
    t_last     = TSYS2TIME(tsys_start);
    
    // sending python #total number of points
    Serial.println(num_iter);
}

//______________________________________________________________________________
//
//
// Perform Input/Ouput on sensors, send data to serial
//
//______________________________________________________________________________

void performIO(const char *what)
{
    const float local_t = TSYS2TIME(TSYS() - tsys_start);
    Serial.print(local_t) ;
    Serial.print("\t");
    Serial.print(what);
    Serial.println("");
}

//______________________________________________________________________________
//
//
// Main Loop
//
//______________________________________________________________________________
void loop()
{
    
    tsys = TSYS() - tsys_start;
    t    = TSYS2TIME(tsys);
    //Serial.println(t);
    if(iter<points_during_rest)
    {
        //______________________________________________________________________
        //
        //
        // resting
        //
        //______________________________________________________________________
        if(t-t_last>=dt_rest)
        {
            // emit a resting point
            performIO("resting");
            // update
            ++iter;
            t_last = t;
        }
        
    }
    else
    {
        
        if(iter<num_iter)
        {
            
            if(resting)
            {
                t_swim  = t;
                resting = false;
            }
            
            //__________________________________________________________________
            //
            //
            // swimming
            //
            //__________________________________________________________________
            if(t-t_last>=dt_swim)
            {
                // initialize t_swim for servo position
               
                // emit a swimming point
                performIO("swimming");
                
                // update
                ++iter;
                t_last = t;
            }
        }
        else
        {
            //__________________________________________________________________
            //
            //
            // Dead...
            //
            //__________________________________________________________________
            
        }
        
    }
    
    
}
