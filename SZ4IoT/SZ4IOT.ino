#include "main.h"
#include "SPIFFS.h"


void SZfunc_float(float *data, int nbEle){
   
      int r1=0,r2=0,r3=0,r4=0,r5=0;
      r1=nbEle;
      size_t outSize;
      unsigned char *bytes;
 
      unsigned long ms1 = micros ();
      bytes = SZ_compress(SZ_FLOAT, data, &outSize, r5, r4, r3, r2, r1);    
      unsigned long ms2 = micros ();
      Serial.print(ms2-ms1);
      Serial.print("\t\t");
      Serial.print(outSize);
      Serial.print("\t");
      size_t byteLength=outSize;
      unsigned long dms1 = micros ();
      float *data2 = (float*) SZ_decompress(SZ_FLOAT, bytes, byteLength, r5, r4, r3, r2, r1);
      unsigned long dms2 = micros ();
      Serial.print(dms2-dms1);
      Serial.print("\t\t");
      
      float *error=(float*)malloc(nbEle*sizeof(float));
      for(int i=0;i<nbEle;i++)
          error[i]=fabs(data[i]-data2[i]);
      float myrms=rms4float(error,nbEle);
      Serial.print(myrms);
      Serial.print("\t");

      float Compratio = ((float)(nbEle*(sizeof(float)))/outSize); 
      Serial.print(Compratio);
      Serial.print("\t");
      Serial.println();
  /*
    for(int i=0;i<10;i++) {
      printf("original %f after com/decom %f\n",data[i],data2[i]);
    }
  */
    free(data2);
    free(bytes);
    free(error);
    
}


void SZfunc_double(double *data, int nbEle){
    
    int r1=0,r2=0,r3=0,r4=0,r5=0;
    r1=nbEle;
    size_t outSize;
    unsigned char *bytes;
 
    unsigned long ms1 = micros ();
    bytes = SZ_compress(SZ_DOUBLE, data, &outSize, r5, r4, r3, r2, r1);    
    unsigned long ms2 = micros ();
    Serial.print(ms2-ms1);
    Serial.print("\t\t");

    Serial.print(outSize);
    Serial.print("\t");
       
    size_t byteLength=outSize;
    unsigned long dms1 = micros ();
    double *data2 = (double*) SZ_decompress(SZ_DOUBLE, bytes, byteLength, r5, r4, r3, r2, r1);
    unsigned long dms2 = micros ();
      
    Serial.print(dms2-dms1);
    Serial.print("\t\t");
      
    double *error=(double*)malloc(nbEle*sizeof(double));
  
    for(int i=0;i<nbEle;i++)
       error[i]= abs(data[i]-data2[i]);
    float myrms=rms4double(error,nbEle);
    Serial.print(myrms);
    Serial.print("\t");

    float Compratio = ((float)(nbEle*(sizeof(double)))/outSize); 
    Serial.print(Compratio);
    Serial.print("\t");
    Serial.println();
  /*
    for(int i=0;i<10;i++) {
      printf("original %f after com/decom %f   error %f\n",data[i],data2[i],error[i]);
    }
  */
    free(data2);
    free(bytes);
    free(error);
}


void SZfunc_INT8(int8_t *data, int nbEle){

    int r1=0,r2=0,r3=0,r4=0,r5=0;
    r1=nbEle;
    size_t outSize;
    unsigned char *bytes;
 
    unsigned long ms1 = micros ();
    bytes = SZ_compress(SZ_INT8, data, &outSize, r5, r4, r3, r2, r1);    
    unsigned long ms2 = micros ();
      
    Serial.print(ms2-ms1);
    Serial.print("\t\t");
    Serial.print(outSize);
    Serial.print("\t");
        
    size_t byteLength=outSize;
    unsigned long dms1 = micros ();
    int8_t *data2 = (int8_t*) SZ_decompress(SZ_INT8, bytes, byteLength, r5, r4, r3, r2, r1);
    unsigned long dms2 = micros ();
      
    Serial.print(dms2-dms1);
    Serial.print("\t\t");
      
    int8_t *error=(int8_t*)malloc(nbEle*sizeof(int8_t));
  
    for(int i=0;i<nbEle;i++)
        error[i]= abs(data[i]-data2[i]);
    float myrms=rms4int8(error,nbEle);
    Serial.print(myrms);
    Serial.print("\t");

    float Compratio = ((float)(nbEle*(sizeof(int8_t)))/outSize); 
    Serial.print(Compratio);
    Serial.print("\t");
    Serial.println();
  
    free(data2);
    free(bytes);
    free(error);
}


void SZfunc_INT16(int16_t *data, int nbEle){

    int r1=0,r2=0,r3=0,r4=0,r5=0;
    r1=nbEle;
    size_t outSize;
    unsigned char *bytes;
 
    unsigned long ms1 = micros ();
    bytes = SZ_compress(SZ_INT16, data, &outSize, r5, r4, r3, r2, r1);    
    unsigned long ms2 = micros ();
      
    Serial.print(ms2-ms1);
    Serial.print("\t\t");
    Serial.print(outSize);
    Serial.print("\t");
        
    size_t byteLength=outSize;
    unsigned long dms1 = micros ();
    int16_t *data2 = (int16_t*) SZ_decompress(SZ_INT16, bytes, byteLength, r5, r4, r3, r2, r1);
    unsigned long dms2 = micros ();
     
    Serial.print(dms2-dms1);
    Serial.print("\t\t");
      
    int16_t *error=(int16_t*)malloc(nbEle*sizeof(int16_t));
  
    for(int i=0;i<nbEle;i++)
        error[i]= abs(data[i]-data2[i]);  
    float myrms=rms4int16(error,nbEle);
    Serial.print(myrms);
    Serial.print("\t");

    float Compratio = ((float)(nbEle*(sizeof(int16_t)))/outSize); 
    Serial.print(Compratio);
    Serial.print("\t");
    Serial.println();
  /*
    for(int i=0;i<10;i++) {
      printf("original %d after com/decom %d\n",data[i],data2[i]);
    }
  */
   
    free(data2);
    free(bytes);
    free(error);
    
}


void SZfunc_INT32(int32_t *data, int nbEle){

    int r1=0,r2=0,r3=0,r4=0,r5=0;
    r1=nbEle;
    size_t outSize;
    unsigned char *bytes;
 
    unsigned long ms1 = micros ();
    bytes = SZ_compress(SZ_INT32, data, &outSize, r5, r4, r3, r2, r1);    
    unsigned long ms2 = micros ();
      
    Serial.print(ms2-ms1);
    Serial.print("\t\t");
    Serial.print(outSize);
    Serial.print("\t");
        
    size_t byteLength=outSize;
    unsigned long dms1 = micros ();
    int32_t *data2 = (int32_t*) SZ_decompress(SZ_INT32, bytes, byteLength, r5, r4, r3, r2, r1);
    unsigned long dms2 = micros ();
      
    Serial.print(dms2-dms1);
    Serial.print("\t\t");
     
    int32_t *error=(int32_t*)malloc(nbEle*sizeof(int32_t));
    for(int i=0;i<nbEle;i++)
        error[i]= abs(data[i]-data2[i]); 
    float myrms=rms4int32(error,nbEle);
    Serial.print(myrms);
    Serial.print("\t");

    float Compratio = ((float)(nbEle*(sizeof(int32_t)))/outSize); 
    Serial.print(Compratio);
    Serial.print("\t");
    Serial.println();
  
    free(data2);
    free(bytes);
    free(error);
 }

void setup() {
  // put your setup code here, to run once:
 Serial.begin(115200); 
 while(!Serial);
 randomSeed(124);
  delay(500);

 int nbEle = 256;
//Serial.printf("nbEle = %d \n", nbEle);
  
//float *data=(float*)malloc(nbEle*sizeof(float));
// double *data=(double*)malloc(nbEle*sizeof(double));
int8_t *data=(int8_t*)malloc(nbEle*sizeof(int8_t));
 //int16_t *data=(int16_t*)malloc(nbEle*sizeof(int16_t));
 // int32_t *data=(int32_t*)malloc(nbEle*sizeof(int32_t));

  Serial.print("Comptime  ");
  Serial.print("CompresdSize\t");
  Serial.print("DeComptime\t");
  Serial.print("RMS\t");
  Serial.println("Compratio");
  
  SZ_Init(NULL);
  errorBoundMode = REL;
  confparams_cpr->errorBoundMode = errorBoundMode;
  confparams_cpr->relBoundRatio= 0.05;

 if (!SPIFFS.begin(true)) {
    Serial.println("An Error has occurred while mounting SPIFFS");
    return;
  } 
  File root = SPIFFS.open("/");
  File file = root.openNextFile();
  
  while(file) {
     for(int i=0;i<nbEle;i++){
      String str = file.readStringUntil('\n');
      //  data[i] = str.toFloat();
       // data[i] = str.toDouble();
        data[i] = str.toInt();
     //  printf("original %f \n", data[i] );

      }
    file.close();
 //  SZfunc_float(data, nbEle);
  //  SZfunc_double(data, nbEle); 
  SZfunc_INT8(data, nbEle);
 // SZfunc_INT16(data, nbEle);
 // SZfunc_INT32(data, nbEle);

    file = root.openNextFile();
}
 free(data);
   
}

void(* resetFunc) (void) = 0;

void loop() {
  // put your main code here, to run repeatedly:

}
