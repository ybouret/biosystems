//////////////////////////////////////////////////////////////////////////////////
//
// Serial
//
//////////////////////////////////////////////////////////////////////////////////
void setup_Serial()
{
  //////////////////////////////////////////////////////////////////////////////////
  //
  // Initialize serial com
  //
  //////////////////////////////////////////////////////////////////////////////////
  Serial.begin(115200);

  // wait until serial port opens for native USB devices
  while (! Serial) {
    delay(1);
  }
}

//////////////////////////////////////////////////////////////////////////////////
//
// LIDARS
//
//////////////////////////////////////////////////////////////////////////////////
#include <Adafruit_VL53L0X.h>

#define MY_LIDAR_ADDR1 0x31
#define MY_LIDAR_ADDR2 0x32
#define MY_LIDAR_CTRL1 9
#define MY_LIDAR_CTRL2 10
#define MY_LIDAR_DELAY 100

Adafruit_VL53L0X lidar1 = Adafruit_VL53L0X();
Adafruit_VL53L0X lidar2 = Adafruit_VL53L0X();

void setup_LIDARS()
{

  Serial.println(F("Initializing LIDARS"));
  Wire.begin();

  //________________________________________________________________________________
  //
  // shutdown every one
  //________________________________________________________________________________
  pinMode(MY_LIDAR_CTRL1, OUTPUT);
  pinMode(MY_LIDAR_CTRL2, OUTPUT);
  digitalWrite(MY_LIDAR_CTRL1, LOW);
  digitalWrite(MY_LIDAR_CTRL2, LOW);
  delay(MY_LIDAR_DELAY);
  digitalWrite(MY_LIDAR_CTRL1, HIGH);
  digitalWrite(MY_LIDAR_CTRL2, HIGH);
  delay(MY_LIDAR_DELAY);

  //________________________________________________________________________________
  //
  // shutdown LIDAR2, setup LIDAR1
  //________________________________________________________________________________
  digitalWrite(MY_LIDAR_CTRL2, LOW);
  delay(MY_LIDAR_DELAY);
  if ( !lidar1.begin(MY_LIDAR_ADDR1) )
  {
    Serial.print(F("LIDAR1 begin failure on A")); Serial.println(MY_LIDAR_CTRL1);
    while (1);
  }

  //________________________________________________________________________________
  //
  // shutdown LIDAR1, wake up LIDAR2, setup LIDAR2,
  //________________________________________________________________________________
  digitalWrite(MY_LIDAR_CTRL2, HIGH);
  delay(MY_LIDAR_DELAY);
  if ( !lidar2.begin(MY_LIDAR_ADDR2) )
  {
    Serial.print(F("LIDAR2 begin failure on A")); Serial.println(MY_LIDAR_CTRL2);
    while (1);
  }

  //________________________________________________________________________________
  //
  // wake up LIDAR1
  //________________________________________________________________________________
  //digitalWrite(MY_LIDAR_CTRL1, HIGH);
  //delay(MY_LIDAR_DELAY);
}


void give_LIDAR(struct Adafruit_VL53L0X *lidar, const int indx)
{
  VL53L0X_RangingMeasurementData_t measure;
  lidar->rangingTest(&measure, false);
  //Serial.print(F("Distance@")); Serial.print(indx,HEX); Serial.print(F(" : "));
  Serial.print(F(" "));
  if (measure.RangeStatus != 4) {  // phase failures have incorrect data
    Serial.print(measure.RangeMilliMeter);
  }
  else
  {
    Serial.print(F("out of range"));
  }
}


void give_LIDARS()
{
  Serial.print(F("Distances: "));
  give_LIDAR(&lidar1, MY_LIDAR_CTRL1);
  give_LIDAR(&lidar2, MY_LIDAR_CTRL2);
  Serial.println(F(""));
}



//////////////////////////////////////////////////////////////////////////////////
//
// 9DOF libraries
//
//////////////////////////////////////////////////////////////////////////////////
#include <Adafruit_Sensor.h>
#include <Adafruit_FXOS8700.h>
Adafruit_FXOS8700 accelmag = Adafruit_FXOS8700(0x8700A, 0x8700B);

void setup_9DOF()
{
  /* Initialise the 9DOF */
  if (!accelmag.begin(ACCEL_RANGE_4G))
  {
    /* There was a problem detecting the FXOS8700 ... check your connections */
    Serial.println("Ooops, no FXOS8700 detected ... Check your wiring!");
    while (1);
  }
}

void give_9DOF()
{
  // 9D0F
  sensors_event_t aevent, mevent;

  /* Get a new sensor event */
  accelmag.getEvent(&aevent, &mevent);


  // Display the accel results (acceleration is measured in m/s^2)
  Serial.print("A ");
  Serial.print("X: "); Serial.print(aevent.acceleration.x, 4); Serial.print("  ");
  Serial.print("Y: "); Serial.print(aevent.acceleration.y, 4); Serial.print("  ");
  Serial.print("Z: "); Serial.print(aevent.acceleration.z, 4); Serial.print("  ");
  Serial.println("m/s^2");

  // Display the mag results (mag data is in uTesla)
  Serial.print("M ");
  Serial.print("X: "); Serial.print(mevent.magnetic.x, 1); Serial.print("  ");
  Serial.print("Y: "); Serial.print(mevent.magnetic.y, 1); Serial.print("  ");
  Serial.print("Z: "); Serial.print(mevent.magnetic.z, 1); Serial.print("  ");
  Serial.println("uT");
}

//////////////////////////////////////////////////////////////////////////////////
//
// Global Setup
//
//////////////////////////////////////////////////////////////////////////////////

void setup()
{
  setup_Serial();
  setup_LIDARS();
  setup_9DOF();
}

void loop()
{
  give_LIDARS();
  give_9DOF();
  delay(200);
}

