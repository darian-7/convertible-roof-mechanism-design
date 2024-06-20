void setup() {
  // Configuration of ports as either inputs or outputs
  pinMode (1, OUTPUT);
  pinMode (3, OUTPUT);
  pinMode (2, INPUT) ;
  pinMode (4, INPUT) ;
  pinMode (5, INPUT) ;
  pinMode (6, INPUT) ;
  pinMode (7, INPUT) ;
  pinMode (8, INPUT) ;
  pinMode (9, INPUT) ;
  pinMode (10, INPUT);
  Serial.begin(9600) ;
}

int ultrasonic = 7; // Initialize pin D7 as input from ultrasonic sensor
int duration = 0; // Sets initial signal duration from ultrasonic sensor to 0s
int distance = (duration/2)/28.5; // Distance formula that uses MEASURED signal duration from ultrasonic sensor
int capacitive = 10; // Initialize pin A2 as input from capacitive sensor 
int capacitance = 0; // Variable that reads the value from the capacitive sensor
int airflow = 9; // Initialize pin A1 as input from the air flow sensor
int halleffect = 8; // Initialize pin D8 as input from hall effect sensor
float revolutions=0; // Sets initial wheel revolutions to 0
int rpm=0; // Sets initial wheel rpm to 0
long  startTime=0; // Sets start time for hall effect sensor tachometer to 0s
long  elapsedTime; // Variable for total elapsed time of hall effect sensor tachometer

void loop() {
  duration = pulseInLong(ultrasonic, HIGH); // Measures time it takes for one pulse from the ultrasonic sensor
  capacitance = analogRead(capacitive); // Reads capacitance value from capacitive sensor
  if (digitalRead(2) == HIGH) { 
    digitalWrite(1, HIGH) ;
    analogWrite(3, 225); } // Turns on motor if switch is pressed to open the roof
  else if (digitalRead(4) == HIGH) {
    digitalWrite(1, LOW) ;
    analogWrite(3, 255) ; } // Turns on motor if switch is pressed to close the roof
  int windADunits = analogRead(airflow); // Reads airflow value from air flow sensor
  float windMPH = pow((((float)windADunits - 264.0) / 85.6814), 3.36814); // Wind formula that reads input from air flow sensor and outputs wind speed in mph
  startTime = millis();
  rpm = (max(1, revolutions) * 60000) / millis()-startTime; // Formula that calculates rpm using hall effect sensor readings and elapsed time of wheel revolutions
  if (windMPH >= 20 || rpm >= 395 || distance <= 25.3 || capacitance >= 1 || digitalRead(5) == HIGH || digitalRead(6) == HIGH) {
    analogWrite(3, 0) ; } // Turns off motor if wind speed >= 20mph OR wheel rotations >= 395 RPM (car speed: 20mph) OR distance above car is <= than 253mm OR capacitance >= 1 microFarad OR if either end stop is reached
}
