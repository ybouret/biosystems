#include <Servo.h>

//______________________________________________________________________________
//
//
// Paramètres modifiables pour les mesures
//
//______________________________________________________________________________

int             Iter_rest          = 5;
int             amplitude          = 100;
float           frequency          = 1;        // 2pi/omega
int             nb_point           = 100;      // number of info points
//int             nb_cycle           = 5;        // swimming #cycles
int nb_cycle=0;
int             Iter_move          = nb_point * nb_cycle; // assign one position per 'move'
float           interFrame_period  = 1000/frequency/nb_point; // time between two points in ms
unsigned int    servo_period       = 15000;                   // microseconds
int             angle_shift        = 0;                       // servo neutral position correction

//______________________________________________________________________________
//
//
// Pin utilisés
//
//______________________________________________________________________________
int      pin_servo          = 9;
int      pin_interFrame     = 12;
int      pin_forceSensor    = 1;
int      pin_tensionEntree  = 3;
int      pin_tensionSortie  = 2;

//______________________________________________________________________________
//
//
// Définition des autres variables
//
//______________________________________________________________________________
unsigned long     baudrate       = 115200;
int               maxIter        = Iter_rest+Iter_move;
int               voltage;
const float       pi             = 3.1415926535;
Servo             myservo;
const int         angle_init     = 90-angle_shift;
unsigned int      last_servo     = 0;
float             sweep          = 2*pi*frequency; // omega
int               angle          = angle_init;
float             ang;
float             dt;
unsigned int      count;
unsigned int      last_interFrame;
bool              flag_interFrame = false; // signal to serial not to record picture
float             tension1;
float             tension2;
float             tension;
unsigned int      t               = micros();
unsigned int      t0              = millis();

bool              init_start=true;
unsigned int      t_start_mu=0;

//______________________________________________________________________________
//
//
// Initialisation
//
//______________________________________________________________________________
void setup()
{
    Serial.begin(baudrate);
    Serial.println(maxIter);
    
    myservo.attach(pin_servo,900,2100);
    myservo.write(angle_init);
    delay(1000);
    
    pinMode(pin_interFrame,OUTPUT);
    digitalWrite(pin_interFrame,LOW);
    
    count           = 0;
    last_interFrame = millis();
    last_servo      = micros();
    t  = millis();
    t0 = millis();
}

//______________________________________________________________________________
//
//
// Serial Sommunication
//
//______________________________________________________________________________
void performIO(const int pinInterFrameValue)
{
    digitalWrite(pin_interFrame,pinInterFrameValue);
    voltage  = analogRead(pin_forceSensor);
    tension1 = analogRead(pin_tensionEntree);
    tension2 = analogRead(pin_tensionSortie);
    tension  = tension2-tension1;
    
    
    Serial.print(micros()-t_start_mu);
    Serial.print("\t");
    Serial.print(voltage);
    Serial.print("\t");
    Serial.print(angle);
    Serial.print("\t");
    Serial.print(flag_interFrame);
    Serial.print("\t");
    Serial.println(tension);
}

//______________________________________________________________________________
//
//
// Main Loop
//
//______________________________________________________________________________
void loop()
{
    if(init_start)
    {
        t_start_mu = micros();
        Serial.println(t_start_mu);
        //Serial.println(1000*interFrame_period*0.5);
        init_start = false;
        last_interFrame = millis();
    }
    
    if(count<=maxIter)
    {
        //______________________________________________________________________
        //
        //
        // Alive
        //
        //______________________________________________________________________
        t = millis();
        if(count>Iter_rest)
        {
            //__________________________________________________________________
            //
            // Swimming
            //__________________________________________________________________
            if(t-last_interFrame>interFrame_period*0.5)
            {
                flag_interFrame =!flag_interFrame;
                if(flag_interFrame)
                {
                    performIO(HIGH);
                    count++;
                }
                else
                {
                    digitalWrite(pin_interFrame,LOW);
                }
                last_interFrame = t;
            }
            if(t*1000-last_servo<servo_period)
            {
                dt    = (millis()-t0)*0.001;
                ang   = 0.5*amplitude*sin(sweep*dt)+90-angle_shift;
                angle = (int) ang;
                myservo.write(angle);
                last_servo = t*1000;
            }
        }
        else
        {
            //__________________________________________________________________
            //
            // Resting: no camera, physical constants out
            //__________________________________________________________________
            if(t-last_interFrame>interFrame_period*0.5)
            {
                flag_interFrame = false; // no picture
                angle = angle_init;
                myservo.write(angle);
                Serial.println(micros()-t_start_mu);
                //performIO(LOW);
                count++;
                last_interFrame = t;
            }
        }
    }
    else
    {
        //______________________________________________________________________
        //
        //
        // Dead
        //
        //______________________________________________________________________
        Serial.end();                // close serial comm
        myservo.write(angle_init);   // go back to sleep
        delay(1000);                 // wait for position to be reached
        myservo.detach();            // done
    }
}
