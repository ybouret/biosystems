#include <Servo.h>         

/* Paramètres modifiables pour les mesures */
int             Iter_rest          = 5;                         
int             amplitude          = 100;     
float           frequency          = 1;     
int             nb_point           = 100;
int             nb_cycle           = 5;
int             Iter_move          = nb_point * nb_cycle;
float           interFrame_period  = 1000/frequency/nb_point;    
unsigned int    servo_period       = 15000;   
int             angle_shift        = 0;

/* Pin utilisés */
int      pin_servo          = 9;
int      pin_interFrame     = 12;
int      pin_forceSensor    = 1;
int      pin_tensionEntree  = 3;
int      pin_tensionSortie  = 2;

/* Définition des autres variables */
unsigned long     baudrate       = 115200;
int               maxIter        = Iter_rest+Iter_move;
int               voltage;
const float       pi             = 3.1415926535;
Servo             myservo;
const int         angle_init     = 90-angle_shift;                 
unsigned int      last_servo     = 0;                  
float             sweep          = 2*pi*frequency;
int               angle          = angle_init;
float             ang;
float             dt;
unsigned int      count;
unsigned int      last_interFrame;
bool              flag_interFrame = false;
float             tension1;
float             tension2; 
float             tension;
unsigned int      t               = micros();
unsigned int      t0              = millis();

/* Initialisation */
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
}

/* Programme principal */
void loop()
{
  if(count<=maxIter)
  {
    t = millis();
    if(count>Iter_rest)
    {
    if(t-last_interFrame>interFrame_period*0.5)
    {
      flag_interFrame =!flag_interFrame;
      if(flag_interFrame)
      {
        digitalWrite(pin_interFrame,HIGH);
        voltage = analogRead(pin_forceSensor);
        
        tension1 = analogRead(pin_tensionEntree);              
        tension2 = analogRead(pin_tensionSortie);             
        tension = tension2-tension1;

        Serial.print(micros());
        Serial.print("\t");
        Serial.print(voltage);
        Serial.print("\t");
        Serial.print(angle);
        Serial.print("\t");
        Serial.print(flag_interFrame);
        Serial.print("\t");
        Serial.println(tension);

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
      dt = (millis()-t0)*0.001;
      ang = 0.5*amplitude*sin(sweep*dt)+90-angle_shift;
      angle = (int) ang;
      myservo.write(angle);
      last_servo = t*1000;
    }
    }
    else
    {
      if(t-last_interFrame>interFrame_period*0.5)
      {
        flag_interFrame =false;
        angle = angle_init;
        
        digitalWrite(pin_interFrame,LOW);
        voltage = analogRead(pin_forceSensor);
        myservo.write(angle);
        tension1 = analogRead(pin_tensionEntree);             
        tension2 = analogRead(pin_tensionSortie);              
        tension = tension2-tension1;

        Serial.print(micros());
        Serial.print("\t");
        Serial.print(voltage);
        Serial.print("\t");
        Serial.print(angle);
        Serial.print("\t");
        Serial.print(flag_interFrame);
        Serial.print("\t");
        Serial.println(tension);

        count++;
        last_interFrame = t;
      }
    }
  }
  else
  {
    Serial.end();
    myservo.write(angle_init);
    delay(1000);
    myservo.detach();
  }
}
