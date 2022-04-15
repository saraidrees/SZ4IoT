# SZ4IoT

SZ4IoT: Error-bounded Lossy Compressor library for IoT devices. 

This library contains five main functions to compress and decompress data. 
* SZ4FLOAT: To compress & decompress the given sample of float data.
* SZ4DOUBLE: To compress & decompress the given sample of double data.
* SZ4INT8: To compress & decompress the given sample of int8_t data.
* SZ4INT16: To compress & decompress the given sample of int16_t data.
* SZ4INT32: To compress & decompress the given sample of int32_t data.

The library was tested and worked well with ESP32, Teensy4.0, and RP2040.

# Installation 

1.	Clone SZ4IoT directory under /Arduino/libraries directory.
2.	Add #include <SZ4IoT.h> to your sketch.
3.	For compressing and decompressing (Float, Double, int8_t, int16_t, and int32_t) data, declare the array and its size and then call it inside the setup () function.
4.	For example, to compress and decompress floating data:

```
#include "SZ4IoT.h"
float data[]={};
void setup() {
  Serial.begin(115200); 
  while(!Serial);
  delay(500);
  int nbEle = sizeof(data) / sizeof(float);
  SZ4FLOAT(data, nbEle);
}
void loop() { }
```
## Credits 
This library is created based on the work done in:

* https://github.com/szcompressor/SZ
