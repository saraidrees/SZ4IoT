
#ifndef SZFORIoT_H
#define SZFORIoT_H

#ifdef __cplusplus
extern "C" {
#endif


void SZ4FLOAT(float *data, int nbEle);
void SZ4DOUBLE(double *data, int nbEle);
void SZ4INT8(int8_t *data, int nbEle);
void SZ4INT16(int16_t *data, int nbEle);
void SZ4INT32(int32_t *data, int nbEle);


#ifdef __cplusplus
}
#endif

#endif
