#include <Arduino.h>

#define USE_LIDARS 1
#define USE_9DOF   0
#define USE_SERVO  0
#define USE_FORCE  0
#define USE_BLE    1
#define USE_SERIAL 1

//////////////////////////////////////////////////////////////////////////////////
//
// BLE
//
//////////////////////////////////////////////////////////////////////////////////
#if 1 == USE_BLE
#include <SPI.h>
#include <Adafruit_BLE.h>
#include <Adafruit_BluefruitLE_SPI.h>
#include <Adafruit_BluefruitLE_UART.h>

#include "BluefruitConfig.h"

// WARNING: this is specific to the Feather Board
Adafruit_BluefruitLE_SPI ble(BLUEFRUIT_SPI_CS, BLUEFRUIT_SPI_IRQ, BLUEFRUIT_SPI_RST);

// TODO: Check the return code of the functions
void setup_BLE()
{
  ble.begin();
  ble.factoryReset();
  ble.echo(false);
  /* Wait for connection */
  while (! ble.isConnected()) {
    delay(500);
  }
  ble.setMode(BLUEFRUIT_MODE_DATA);
}

void give_BLE()
{
  /*
  const float t = 1.0e-6f * float( micros() );
  ble.print( 90.0f + 90.0f * sin( 6.14f * (t / 10.0f) ) );
  ble.println();
  */
}
#endif

//////////////////////////////////////////////////////////////////////////////////
//
// Serial
//
//////////////////////////////////////////////////////////////////////////////////
#if 1 == USE_SERIAL
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
#endif

#if 1 == USE_LIDARS
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

Adafruit_VL53L0X *lidar1 = NULL; //Adafruit_VL53L0X();
Adafruit_VL53L0X *lidar2 = NULL; //Adafruit_VL53L0X();


void setup_LIDARS()
{


  Serial.print(F("Initializing LIDARS/2*")); Serial.println(sizeof(Adafruit_VL53L0X));
  lidar1 = new Adafruit_VL53L0X();
  lidar2 = new Adafruit_VL53L0X();

  //________________________________________________________________________________
  //
  // shutdown every one
  //________________________________________________________________________________
  pinMode(MY_LIDAR_CTRL1, OUTPUT);
  pinMode(MY_LIDAR_CTRL2, OUTPUT);
  digitalWrite(MY_LIDAR_CTRL1, LOW);
  digitalWrite(MY_LIDAR_CTRL2, LOW);
  delay(MY_LIDAR_DELAY);

  //________________________________________________________________________________
  //
  // wake up everyone: this is a reset...
  //________________________________________________________________________________
  digitalWrite(MY_LIDAR_CTRL1, HIGH);
  digitalWrite(MY_LIDAR_CTRL2, HIGH);
  delay(MY_LIDAR_DELAY);

  //________________________________________________________________________________
  //
  // shutdown LIDAR2, setup LIDAR1
  //________________________________________________________________________________
  digitalWrite(MY_LIDAR_CTRL2, LOW);
  delay(MY_LIDAR_DELAY);
  if ( !lidar1->begin(MY_LIDAR_ADDR1) )
  {
    Serial.println(F("LIDAR1 failure"));
    while (1);
  }

  //________________________________________________________________________________
  //
  // shutdown LIDAR1, wake up LIDAR2, setup LIDAR2,
  //________________________________________________________________________________
  digitalWrite(MY_LIDAR_CTRL2, HIGH);
  delay(MY_LIDAR_DELAY);
  if ( !lidar2->begin(MY_LIDAR_ADDR2) )
  {
    Serial.println(F("LIDAR2 failure"));
    while (1);
  }

  //________________________________________________________________________________
  //
  // wake up LIDAR1
  //________________________________________________________________________________
  //digitalWrite(MY_LIDAR_CTRL1, HIGH);
  //delay(MY_LIDAR_DELAY);
}


void give_LIDAR(struct Adafruit_VL53L0X *lidar)
{
  VL53L0X_RangingMeasurementData_t measure;
  lidar->rangingTest(&measure, false);
#if 1 == USE_BLE
  if (measure.RangeStatus != 4) {  // phase failures have incorrect data
    const int rm = int(measure.RangeMilliMeter);
    ble.print(rm);
  }
  else
  {
    ble.print(F("-1"));
  }
#endif

  Serial.print(F(" "));
  if (measure.RangeStatus != 4) {  // phase failures have incorrect data
    Serial.print(measure.RangeMilliMeter);
  }
  else
  {
    Serial.print(F("out"));
  }
}


void give_LIDARS()
{
  Serial.print(F("Distances    : "));
  give_LIDAR(lidar1);
#if 1 == USE_BLE
  ble.print(F(","));
#endif
  give_LIDAR(lidar2);
#if 1 == USE_BLE
  ble.println();
#endif

  Serial.println();
}

#endif //! USE_LIDARS

//////////////////////////////////////////////////////////////////////////////////
//
// 9DOF libraries
//
//////////////////////////////////////////////////////////////////////////////////
#if 1 == USE_9DOF
#include <Adafruit_Sensor.h>
#include <Adafruit_FXOS8700.h>
Adafruit_FXOS8700 accelmag = Adafruit_FXOS8700(0x8700A, 0x8700B);

void setup_9DOF()
{
  /* Initialise the 9DOF */
  if (!accelmag.begin(ACCEL_RANGE_4G))
  {
    /* There was a problem detecting the FXOS8700 ... check your connections */
    Serial.println(F("Ooops, no FXOS8700 detected ... Check your wiring!"));
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
  Serial.print(F("Accelerations: "));

  Serial.print(F("[ "));
  Serial.print(aevent.acceleration.x, 4); Serial.print(F(" "));
  Serial.print(aevent.acceleration.y, 4); Serial.print(F(" "));
  Serial.print(aevent.acceleration.z, 4);
  Serial.println(F("] m/s^2"));

  // Display the mag results (mag data is in uTesla)
  Serial.print(F("Magnetism    : "));
  Serial.print(F("[ "));
  Serial.print(mevent.magnetic.x, 1); Serial.print(F(" "));
  Serial.print(mevent.magnetic.y, 1); Serial.print(F(" "));
  Serial.print(mevent.magnetic.z, 1);
  Serial.println(F("] uT"));
}
#endif // USE_9DOF

//////////////////////////////////////////////////////////////////////////////////
//
// Force Sensors
//
//////////////////////////////////////////////////////////////////////////////////
#if 1 == USE_FORCE
#define PRESSURE_PIN 0
#define FLEXION_PIN  5

void setup_Forces()
{
  pinMode(PRESSURE_PIN, INPUT);
  pinMode(FLEXION_PIN, INPUT);
}

void give_Forces()
{
  const int pressure = analogRead(PRESSURE_PIN);
  const int flexion  = analogRead(FLEXION_PIN);
  Serial.print(F("Analog Forces:"));
  Serial.print(F(" Pressure=")); Serial.print(pressure);
  Serial.print(F(" Flexion="));  Serial.print(flexion);
  Serial.println(F(""));
}
#endif // USE_FORCE

//////////////////////////////////////////////////////////////////////////////////
//
// Servo, if any
//
//////////////////////////////////////////////////////////////////////////////////
#if 1 == USE_SERVO
#include <Servo.h>
#define SERVO_PIN 13
Servo servo;
void setup_Servo()
{
  servo.attach(SERVO_PIN);
  servo.write(90.0f);
}

void give_Servo()
{
  const float t = 1.0e-6f * float( micros() );
  servo.write( 90.0f + 90.0f * sin( 6.14f * (t / 10.0f) ) );
}
#endif



//////////////////////////////////////////////////////////////////////////////////
//
// Global Setup
//
//////////////////////////////////////////////////////////////////////////////////

void setup()
{
#if 1 == USE_SERIAL
  setup_Serial();
#endif

#if 1 == USE_LIDARS
  setup_LIDARS();
#endif

#if 1 == USE_9DOF
  setup_9DOF();
#endif

#if 1 == USE_FORCE
  setup_Forces();
#endif

#if 1 == USE_SERVO
  setup_Servo();
#endif

#if 1 == USE_BLE
  setup_BLE();
#endif

}

void loop()
{
#if 1 == USE_LIDARS
  give_LIDARS();
#endif

#if 1 == USE_9DOF
  give_9DOF();
#endif

#if 1 == USE_FORCE
  give_Forces();
#endif

#if 1 == USE_SERVO
  give_Servo();
#endif

#if 1 == USE_BLE
  give_BLE();
#endif

  delay(100);
}

