#include <Servo.h>
//______________________________________________________________________________
//
//
// communication variables
//
//______________________________________________________________________________
unsigned long     baudrate       = 115200;


//______________________________________________________________________________
//
//
// Simulation parameters, user defined
//
//______________________________________________________________________________

float    frequency          = 1.0;     // servo frequency, in Hz
float    amplitude          = 100.0;   // in degrees
unsigned num_periods        = 5;       // number of servo periods during motion
unsigned points_per_period  = 20;      // number of saved points per cycle
float    resting_time       = 0.5;     // resting initial time
unsigned points_during_rest = 10;      // points to save during resting time

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
const float       pi  = 3.1415926535;
const float omega     = 2.0*pi*frequency;
float simulation_time = resting_time + num_periods * period;




unsigned long     tsys_start     = 0;      // in microseconds
unsigned long     tsys           = 0;
float             t              = 0;
float             t_last         = 0;
float             t_save         = 0.25;

#define TSYS()       (micros())
#define TSYS2TIME(T) (1e-6* ((float)(T)))

//______________________________________________________________________________
//
//
// Initialisation
//
//______________________________________________________________________________
void setup()
{
    Serial.begin(baudrate);
    
    // initialize
    tsys_start = TSYS();
    t_last     = TSYS2TIME(tsys_start);
    //Serial.println(tsys_start);
    
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
            //__________________________________________________________________
            //
            //
            // swimming
            //
            //__________________________________________________________________
            if(t-t_last>=dt_swim)
            {
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
