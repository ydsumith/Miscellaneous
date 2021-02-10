#include <Wire.h>
#include <LiquidCrystal_I2C.h>
#include "DHT.h"

#define DHTPIN 13 
// Uncomment whatever type you're using!
#define DHTTYPE DHT11   // DHT 11
//#define DHTTYPE DHT22   // DHT 22  (AM2302), AM2321
//#define DHTTYPE DHT21   // DHT 21 (AM2301)
DHT dht(DHTPIN, DHTTYPE);

LiquidCrystal_I2C lcd (0x27,16,2); // set the LCD address to 0x27 for a16 chars and 2 line display

int temperaturepin = 0;
int ledwhite = 2;//8
int ledyellow = 8;
int counter =0;
int Tcur;
int Tset = 70;
int minusbutton = 5;
int plusbutton = 4;

void setup() {
  pinMode (temperaturepin, INPUT); 
  pinMode(ledwhite, OUTPUT);
  pinMode(ledyellow, OUTPUT);
  lcd.init (); // initialize the lcd
  lcd.backlight ();
  lcd.clear (); // Clear the screen
  Serial.begin (9600);
  dht.begin();
}

// the loop function runs over and over again forever
void loop() {
  digitalWrite(ledwhite, LOW);   // turn the LED on (HIGH is the voltage level)
  digitalWrite(ledyellow, HIGH);   // turn the LED on (HIGH is the voltage level)
  
  // Reading temperature or humidity takes about 250 milliseconds!
  // Sensor readings may also be up to 2 seconds 'old' (its a very slow sensor)
  float h = dht.readHumidity();
  // Read temperature as Celsius (the default)
  //float t = dht.readTemperature();
  // Read temperature as Fahrenheit (isFahrenheit = true)
  float f = dht.readTemperature(true);
  Tcur = (int) f;
  
  // Check if any reads failed and exit early (to try again).
  if (isnan(h) || isnan(f)) {
    Serial.println(F("Failed to read from DHT sensor!"));
    return;
  }
  
  //lcd.clear ();
  lcd.setCursor (0,0);
  lcd.print("Room:");
  lcd.setCursor (5,0); // (space, row)
  lcd.print(Tcur);
  lcd.setCursor (7,0);
  lcd.print("F Hum:");
  lcd.setCursor (13,0);
  lcd.print(h);
  lcd.setCursor (15,0);
  lcd.print("%");

  
  lcd.setCursor (0,1);
  lcd.print("Set Temp:"); // LED print keyestudio!
  if(digitalRead(minusbutton)==0){
    Tset = Tset -1;
  }
  if(digitalRead(plusbutton)==0){
    Tset = Tset +1;
  }
  lcd.setCursor (9,1);
  lcd.print(Tset);
  
  delay(500);                       // wait for a second
  digitalWrite(ledwhite, HIGH);    // turn the LED off by making the voltage LOW
  digitalWrite(ledyellow, LOW);   // turn the LED on (HIGH is the voltage level)
  counter = counter +1;
  if (counter == 10){
    counter =0;
    lcd.clear (); // Clear the screen
  }
  delay(500);                       // wait for a second
}
