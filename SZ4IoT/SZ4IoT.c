#include <string.h>
#include "stdint.h"
#include"SZ_lib.h"


sz_params *confparams_cpr;
int  errorBoundMode;

void SZ4FLOAT(float *data, int nbEle){
   
    int r1=0,r2=0,r3=0,r4=0,r5=0;
    r1=nbEle;
    size_t outSize;
    unsigned char *bytes;
    SZ_Init(NULL);
    errorBoundMode = REL;
    confparams_cpr->errorBoundMode = errorBoundMode;
    unsigned long ms1 = micros ();
    bytes = SZ_compress(SZ_FLOAT, data, &outSize, r5, r4, r3, r2, r1);    
    unsigned long ms2 = micros ();
    printf("Time for compressing %d bytes :%d\n",nbEle*sizeof(float),ms2-ms1); 
    printf("COMPRESSED SIZE %d\n",outSize);
      
    size_t byteLength=outSize;
    unsigned long dms1 = micros ();
    float *data2 = (float*) SZ_decompress(SZ_FLOAT, bytes, byteLength, r5, r4, r3, r2, r1);
    unsigned long dms2 = micros ();
    printf("Time for decompressing %d bytes :%d\n",nbEle*sizeof(float),dms2-dms1);
     
    float *error=(float*)malloc(nbEle*sizeof(float));
    for(int i=0;i<nbEle;i++)
        error[i]=fabs(data[i]-data2[i]);
    float myrms=rms4float(error,nbEle);
    printf("RMS %f\n",myrms);
    

    free(data2);
    free(bytes);
    free(error);
    
}


void SZ4DOUBLE(double *data, int nbEle){
    
    int r1=0,r2=0,r3=0,r4=0,r5=0;
    r1=nbEle;
    size_t outSize;
    unsigned char *bytes;
    SZ_Init(NULL);
    errorBoundMode = REL;
    confparams_cpr->errorBoundMode = errorBoundMode;
    unsigned long ms1 = micros ();
    bytes = SZ_compress(SZ_DOUBLE, data, &outSize, r5, r4, r3, r2, r1);    
    unsigned long ms2 = micros ();
      
    printf("Time for compressing %d bytes :%d\n",nbEle*sizeof(double),ms2-ms1); 

    printf("COMPRESSED SIZE %d\n",outSize);
    size_t byteLength=outSize;
    unsigned long dms1 = micros ();
    double *data2 = (double*) SZ_decompress(SZ_DOUBLE, bytes, byteLength, r5, r4, r3, r2, r1);
    unsigned long dms2 = micros ();
      
    printf("Time for decompressing %d bytes :%d\n",nbEle*sizeof(double),dms2-dms1);
      
    double *error=(double*)malloc(nbEle*sizeof(double));
  
    for(int i=0;i<nbEle;i++)
       error[i]= fabs(data[i]-data2[i]);
    float myrms=rms4double(error,nbEle);
    printf("RMS %f\n",myrms);
   
    free(data2);
    free(bytes);
    free(error);
}


void SZ4INT8(int8_t *data, int nbEle){
    
   int r1=0,r2=0,r3=0,r4=0,r5=0;
   r1=nbEle;
   size_t outSize;
   unsigned char *bytes;
   SZ_Init(NULL);
   errorBoundMode = REL;
   confparams_cpr->errorBoundMode = errorBoundMode;
   unsigned long ms1 = micros ();
   bytes = SZ_compress(SZ_INT8, data, &outSize, r5, r4, r3, r2, r1);    
   unsigned long ms2 = micros ();
      
   printf("Time for compressing %d bytes :%d\n",nbEle*sizeof(int8_t),ms2-ms1); 
   printf("COMPRESSED SIZE %d\n",outSize);
    
   size_t byteLength=outSize;
   unsigned long dms1 = micros ();
   int8_t *data2 = (int8_t*) SZ_decompress(SZ_INT8, bytes, byteLength, r5, r4, r3, r2, r1);
   unsigned long dms2 = micros ();
      
   printf("Time for decompressing %d bytes :%d\n",nbEle*sizeof(int8_t),dms2-dms1);
      
   int8_t *error=(int8_t*)malloc(nbEle*sizeof(int8_t));
  
   for(int i=0;i<nbEle;i++)
       error[i]= abs(data[i]-data2[i]);
   float myrms=rms4int8(error,nbEle);
   printf("RMS %f\n",myrms);
    

   free(data2);
   free(bytes);
   free(error);
}


void SZ4INT16(int16_t *data, int nbEle){
  
   int r1=0,r2=0,r3=0,r4=0,r5=0;
   r1=nbEle;
   size_t outSize;
   unsigned char *bytes;
   SZ_Init(NULL);
   errorBoundMode = REL;
   confparams_cpr->errorBoundMode = errorBoundMode;
   unsigned long ms1 = micros ();
   bytes = SZ_compress(SZ_INT16, data, &outSize, r5, r4, r3, r2, r1);    
   unsigned long ms2 = micros ();
      
   printf("Time for compressing %d bytes :%d\n",nbEle*sizeof(int16_t),ms2-ms1); 

   printf("COMPRESSED SIZE %d\n",outSize);
    
   size_t byteLength=outSize;
   unsigned long dms1 = micros ();
   int16_t *data2 = (int16_t*) SZ_decompress(SZ_INT16, bytes, byteLength, r5, r4, r3, r2, r1);
   unsigned long dms2 = micros ();
     
   printf("Time for decompressing %d bytes :%d\n",nbEle*sizeof(int16_t),dms2-dms1);
          
   int16_t *error=(int16_t*)malloc(nbEle*sizeof(int16_t));
  
   for(int i=0;i<nbEle;i++)
        error[i]= abs(data[i]-data2[i]);  
   float myrms=rms4int16(error,nbEle);
   printf("RMS %f\n",myrms);
  
   free(data2);
   free(bytes);
   free(error);
    
}


void SZ4INT32(int32_t *data, int nbEle){

   int r1=0,r2=0,r3=0,r4=0,r5=0;
   r1=nbEle;
   size_t outSize;
   unsigned char *bytes;
   SZ_Init(NULL);
   errorBoundMode = REL;
   confparams_cpr->errorBoundMode = errorBoundMode;
   unsigned long ms1 = micros ();
   bytes = SZ_compress(SZ_INT32, data, &outSize, r5, r4, r3, r2, r1);    
   unsigned long ms2 = micros ();   
   printf("Time for compressing %d bytes :%d\n",nbEle*sizeof(int32_t),ms2-ms1); 
   printf("COMPRESSED SIZE %d\n",outSize);
     
   size_t byteLength=outSize;
   unsigned long dms1 = micros ();
   int32_t *data2 = (int32_t*) SZ_decompress(SZ_INT32, bytes, byteLength, r5, r4, r3, r2, r1);
   unsigned long dms2 = micros ();
      
   printf("Time for decompressing %d bytes :%d\n",nbEle*sizeof(int32_t),dms2-dms1);
       
   int32_t *error=(int32_t*)malloc(nbEle*sizeof(int32_t));
   for(int i=0;i<nbEle;i++)
   error[i]= abs(data[i]-data2[i]); 
   float myrms=rms4int32(error,nbEle);
   printf("RMS %f\n",myrms);

    free(data2);
    free(bytes);
    free(error);

 }
