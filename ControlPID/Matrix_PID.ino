/*
 *  CONTROL PID MATRIZ DE LED - © Mecatronium Chips 2020-2021 - Victoriano Montemayor 
 *  victoriano@mecatronium.com
 *  v2.0 - Ene 2021 - Simplificacion pines
 *  
 *  Este programa controla la intensidad de luz medida por la Fotoresistencia por medio 
 *  de un control PID, con calores de Kp, Ki, Kd configurables
 *  Para el kit de Control de Iluminación con Fotoresistencia y Matriz de LED 8x8
 *  
 *  CONEXIONES
 *  
 *  pin A0 se conecta a Potenciometro
 *  pin A1 se conecta a Fotoresistencia
 *  
 *  pin 10 se conecta a CLK
 *  pin 11 se conecta a CS 
 *  pin 12 se conecta a DIN
 * 
 */


#include "LedControl.h"
#include <PID_v1.h>

#define potPin A1         //Pin al que se conecta el Potenciómetro
#define sensorPin A0      //Pin al que se conecta la Fotoresistencia

#define NUMVALORES 100    //Numero de valores a tomar para promediar
#define MS_DELAY 1000     //Tiempo en ms a esperar antes de empezar control (tiempo para estabilizar medicion de sensor)
#define MAX_SETPOINT 100  //Valor máximo permitido de Setpoint (cuando el potenciómetro este al 100%)

int sensorValue = 0;    //Variable para guardar la lectura de la Fotoresistencia
int values[NUMVALORES]; //Arreglo con los ultimos NUMVALORES valores de la Fotoresistencia
int potValues[10];      //Arreglo con los ultimos 10 valores del Potenciómetro
double Setpoint, Input, Output;  //Variables para control PID

//Valores de k's para PID
double Kp = 3.278;    //Proporcional
double Ki = 0.120;    //Integral
double Kd = 0;    //Derivativa

LedControl lc = LedControl(12,10,11,1); //Inicializar Control de Matriz de LED
PID myPID(&Input, &Output, &Setpoint, Kp, Ki, Kd, DIRECT);  //Inicializar control PID. Documentación: https://playground.arduino.cc/Code/PIDLibraryConstructor/


//Inicializar arreglo en 0's
void initValues(){
  for(int i = 0; i < NUMVALORES; i++){
    values[i] = 0;
  }
}

int readPot(int potPin);

void setup() {

  pinMode(A0, INPUT); //Configurar como entrada
  pinMode(A1, INPUT); //Configurar como entrada
  Serial.begin(9600); //Inicializar puerto serial
  
  lc.shutdown(0,false); //Despertar Matriz LED
  lc.setIntensity(0,15); //Establecer intensidad de LED
  lc.clearDisplay(0);   //Apagar todos los LED

  initValues();   //Todo a 0

  myPID.SetMode(AUTOMATIC);               //Modo de control Automático
  myPID.SetOutputLimits(0,64);            //Los valores mínimos y máximos para el Output (0 a 64 LEDs)
  delay(MS_DELAY);                        //Esperar a que se estabilice sensor
}

void loop() {

  int valorPot = readPot(potPin); //Leer valor del potenciómetro
  double valorSensor = readValue();     //Leer valor del sensor
     
  Setpoint = calcularSetpoint(valorPot);  //Calcular Setpoint
  Input = ((valorSensor)/820*100) ;   //Establecer como Input como porcentaje relativo max
  myPID.Compute();                        //Calcular valor de Control > Se actualiza variable Output
  prendeLeds(Output);                     //Enviar valor de Output (Prender LEDs correspondientes)[0 - 64]

  //Graficar
  Serial.print(Setpoint);
  Serial.print(",");
  Serial.print(Output/64*100);
  Serial.print(",");
  Serial.print(Input);
  Serial.print(",");
  Serial.print(0);
  Serial.print(" ");
  Serial.print(",");
  Serial.print(100);
  Serial.println(" ");
  
  delay(1);
}

//Leer el valor actual del sensor y promediar con los ultimos NUMVALORES valores
//para crear un valor mas estable
int readValue(){
  //Recorrer valores
   for(int i = 0; i < NUMVALORES - 1; i++){
    values[i] = values[i+1];
  }

  //Leer el mas nuevo
  values[NUMVALORES - 1] = analogRead(sensorPin);
  
  //Calcular y regresar promedio
  long suma = 0;
  for(int i = 0; i < NUMVALORES; i++){
    suma += values[i];
  }
  double promedio = suma / NUMVALORES;
  
  return promedio;
}


int readPot(int potPin){
  //Recorrer valores
  int sum = 0;
   for(int i = 1; i < 100; i++){
    sum = sum + analogRead(potPin);
  }

  double mean = sum / 100;
  
  return mean;
}

//Prender el numeroDeLeds indicado
void prendeLeds(int numeroDeLeds){
  for(int i = 0; i < 8; i++){
    for (int j = 0; j < 8; j++){
      if( (i*8) + j <= numeroDeLeds){
        lc.setLed(0, i, j, true); 
      } else {
        lc.setLed(0, i, j, false);
      }
    }
  }
}

//Convertir 0% a 100%  >>  0 a 64 LED
int calcularLEDS(int valorEnPorcentaje){
  return (valorEnPorcentaje * 64) / 100;
}

//Convertir valor de Pot (0 a 1023) >> Valor Setpoint (0 a MAX_SETPOINT)
double calcularSetpoint(double valorPotenciometro){
  return (valorPotenciometro * MAX_SETPOINT * 1.0) / 1023;
}
