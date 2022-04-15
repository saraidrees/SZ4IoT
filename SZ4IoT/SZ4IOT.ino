#include <stdio.h> 
#include "SZ4IOT.h"

int32_t data[]={};

void setup() {
  // put your setup code here, to run once:
  Serial.begin(115200); 
  while(!Serial);
  delay(500);

  int nbEle = sizeof(data) / sizeof(int32_t);
  Serial.printf("nbEle = %d\n",nbEle);
    
// SZ4FLOAT(data, nbEle);
// SZ4DOUBLE(data, nbEle);
 //SZ4INT8(data, nbEle);
 //SZ4INT16(data, nbEle);
 SZ4INT32(data, nbEle);
 
}
void loop() {
  // put your main code here, to run repeatedly:
}
