
#ifndef SZ_LIB_H
#define SZ_LIB_H

#include "stdint.h"
#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif



#define LITTLE_ENDIAN_DATA 0 //refers to the endian type of the data read from the disk
#define BIG_ENDIAN_DATA 1 //big_endian (ppc, max, etc.) ; little_endian (x86, x64, etc.)

#define LITTLE_ENDIAN_SYSTEM 0 //refers to the endian type of the system
#define BIG_ENDIAN_SYSTEM 1
#define PASTRI 103
#define HZ 102
#define SZ 101

#define SZ_NO_REGRESSION 0   
#define SZ_WITH_LINEAR_REGRESSION 1

#define MIN_NUM_OF_ELEMENTS 20
#define MIN_ZLIB_DEC_ALLOMEM_BYTES 100000
#define SZ_ZLIB_BUFFER_SIZE 1024

#define SZ_BEST_SPEED 0
#define SZ_BEST_COMPRESSION 1
#define SZ_DEFAULT_COMPRESSION 2
#define SZ_TEMPORAL_COMPRESSION 3

#define GZIP_COMPRESSOR 0 //i.e., ZLIB_COMPRSSOR
#define ZSTD_COMPRESSOR 1

#define PW_REL 10

#define SZ_FLOAT 0
#define SZ_INT8 1
#define SZ_INT16 2
#define SZ_INT32 3
#define SZ_DOUBLE 4


#define ABS 0
#define REL 1
#define ABS_AND_REL 2
#define ABS_OR_REL 3
#define PSNR 4

#define SZ_PWR_MIN_TYPE 0
#define SZ_PWR_AVG_TYPE 1
#define SZ_PWR_MAX_TYPE 2

//SUCCESS returning status
#define SZ_SCES 0  //successful
#define SZ_NSCS -1 //Not successful
#define SZ_FERR -2 //Failed to open input file
#define SZ_TERR -3 //wrong data type (should be only float or double)
#define SZ_DERR -4 //dimension error
#define SZ_MERR -5 //sz_mode error
#define SZ_BERR -6 //bound-mode error (should be only ABS, REL, ABS_AND_REL, ABS_OR_REL, or PW_REL)

#define NORM 5

#define SZ_INT8_MIN -128
#define SZ_INT8_MAX 127
#define SZ_INT16_MIN -32768
#define SZ_INT16_MAX 32767
#define SZ_INT32_MIN -2147483648
#define SZ_INT32_MAX 2147483647
#define MetaDataByteLength_double 36 //meta data length for double type
#define SZ_Transpose 104
#define DynArrayInitLen 1024


#define computeMinMax(data) \
        for(i=1;i<size;i++)\
        {\
                data_ = data[i];\
                if(min>data_)\
                        min = data_;\
                else if(max<data_)\
                        max = data_;\
        }\


typedef struct sz_params
{
  int dataType;
  unsigned int max_quant_intervals; //max number of quantization intervals for quantization
  unsigned int quantization_intervals; 
  unsigned int maxRangeRadius;
  int sol_ID;// it's always SZ, unless the setting is PASTRI compression mode (./configure --enable-pastri)
  int losslessCompressor;
  int sampleDistance; //2 bytes
  float predThreshold;  // 2 bytes
  int szMode; //* 0 (best speed) or 1 (better compression with Gzip) or 3 temporal-dimension based compression
  int gzipMode; //* four options: Z_NO_COMPRESSION, or Z_BEST_SPEED, Z_BEST_COMPRESSION, Z_DEFAULT_COMPRESSION
  int  errorBoundMode; //4bits (0.5byte), //ABS, REL, ABS_AND_REL, or ABS_OR_REL, PSNR, or PW_REL, PSNR
  double absErrBound; //absolute error bound
  double relBoundRatio; //value range based relative error bound ratio
  double psnr; //PSNR
  double pw_relBoundRatio; //point-wise relative error bound
  int segment_size; //only used for 2D/3D data compression with pw_relBoundRatio
  int pwr_type; //only used for 2D/3D data compression with pw_relBoundRatio
  
  int snapshotCmprStep; //perform single-snapshot-based compression if time_step == snapshotCmprStep
  int predictionMode;
    int withRegression;

  int randomAccess;
  int protectValueRange; //0 or 1

  int accelerate_pw_rel_compression;

  
  double dmin, dmax;
  double normErr;


} sz_params;

typedef struct node_t {
  struct node_t *left, *right;
  size_t freq;
  char t; //in_node:0; otherwise:1
  unsigned int c;
} *node;

typedef struct HuffmanTree {
  int stateNum;
  int allNodes;
  struct node_t* pool;
  node *qqq, *qq; //the root node of the HuffmanTree is qq[1]
  int n_nodes; //n_nodes is for compression
  int qend; 
  uint64_t **code;
  unsigned char *cout;
  int n_inode; //n_inode is for decompression
} HuffmanTree;


typedef struct sz_exedata
{
  char optQuantMode;  //opt Quantization (0: fixed ; 1: optimized)  
  int intvCapacity; // the number of intervals for the linear-scaling quantization
  int intvRadius;  // the number of intervals for the radius of the quantization range (intvRadius=intvCapacity/2)
  int SZ_SIZE_TYPE; //the length (# bytes) of the size_t in the system at runtime //4 or 8: sizeof(size_t) 
} sz_exedata;

typedef union lfloat
{
    float value;
    unsigned int ivalue;
    unsigned char byte[4];
} lfloat;
typedef union ldouble
{
    double value;
    long long lvalue;
    unsigned char byte[8];
} ldouble;

typedef struct DynamicIntArray
{ 
  unsigned char* array; //char* (one byte) is enough, don't have to be int*
  size_t size;
  size_t capacity;
} DynamicIntArray;


 
typedef struct TightDataPointStorageF
{
  size_t dataSeriesLength;
  int allSameData;
  double realPrecision; //it's used as the pwrErrBoundRatio when errBoundMode==PW_REL
  float medianValue;
  char reqLength;
  char radExpo; //used to compute reqLength based on segmented precisions in "pw_rel_compression"
  
  int stateNum;
  int allNodes;
  
  size_t exactDataNum;
  float reservedValue;
  
  unsigned char* rtypeArray;
  size_t rtypeArray_size;
  
  float minLogValue;

  unsigned char* typeArray; //its size is dataSeriesLength/4 (or xxx/4+1) 
  size_t typeArray_size;
  
  unsigned char* leadNumArray; //its size is exactDataNum/4 (or exactDataNum/4+1)
  size_t leadNumArray_size;
  
  unsigned char* exactMidBytes;
  size_t exactMidBytes_size;
  
  unsigned char* residualMidBits;
  size_t residualMidBits_size;
  
  unsigned int intervals; //quantization_intervals
  
  unsigned char isLossless; //a mark to denote whether it's lossless compression (1 is yes, 0 is no)
  
  size_t segment_size;
  
  unsigned char* pwrErrBoundBytes;
  int pwrErrBoundBytes_size;
  
  unsigned char* raBytes;
  size_t raBytes_size;
  
} TightDataPointStorageF;


typedef struct TightDataPointStorageI
{
  size_t dataSeriesLength;
  int allSameData;
  double realPrecision; //it's used as the pwrErrBoundRatio when errBoundMode==PW_REL
  size_t exactDataNum;
  long minValue;
  int exactByteSize;
  int dataTypeSize; //the size of data type, e.g., it's 4 when data type is int32_t
  
  int stateNum;
  int allNodes;
  
  unsigned char* typeArray; //its size is dataSeriesLength/4 (or xxx/4+1) 
  size_t typeArray_size;
  
  unsigned char* exactDataBytes;
  size_t exactDataBytes_size;
  
  unsigned int intervals; //quantization_intervals
  
  unsigned char isLossless; //a mark to denote whether it's lossless compression (1 is yes, 0 is no)

} TightDataPointStorageI;

typedef struct TightDataPointStorageD
{
  size_t dataSeriesLength;
  int allSameData;
  double realPrecision;
  double medianValue;
  char reqLength; 
  char radExpo; //used to compute reqLength based on segmented precisions in "pw_rel_compression"

  double minLogValue;

  int stateNum;
  int allNodes;

  size_t exactDataNum;
  double reservedValue;
  
  unsigned char* rtypeArray;
  size_t rtypeArray_size;
  
  unsigned char* typeArray; //its size is dataSeriesLength/4 (or xxx/4+1) 
  size_t typeArray_size;
  
  unsigned char* leadNumArray; //its size is exactDataNum/4 (or exactDataNum/4+1)
  size_t leadNumArray_size;
  
  unsigned char* exactMidBytes;
  size_t exactMidBytes_size;
  
  unsigned char* residualMidBits;
  size_t residualMidBits_size;
  
  unsigned int intervals;
  
  unsigned char isLossless; //a mark to denote whether it's lossless compression (1 is yes, 0 is no)
  
  size_t segment_size;
  
  unsigned char* pwrErrBoundBytes;
  int pwrErrBoundBytes_size;
    
  unsigned char* raBytes;
  size_t raBytes_size;
  
  unsigned char plus_bits;
  unsigned char max_bits;
  
} TightDataPointStorageD;

typedef struct DynamicByteArray
{ 
  unsigned char* array;
  size_t size;
  size_t capacity;
} DynamicByteArray;

typedef struct FloatValueCompressElement
{
   float data;
  int curValue;
  unsigned char curBytes[4]; //big_endian
  int reqBytesLength;
  int resiBitsLength;
} FloatValueCompressElement;

typedef struct DoubleValueCompressElement
{
  double data;
  long curValue;
  unsigned char curBytes[8]; //big_endian
  int reqBytesLength;
  int resiBitsLength;
} DoubleValueCompressElement;

typedef struct LossyCompressionElement
{
  int leadingZeroBytes; //0,1,2,or 3
  unsigned char integerMidBytes[8];
  int integerMidBytes_Length; //they are mid_bits actually
  //char curBytes[8];
  //int curBytes_Length; //4 for single_precision or 8 for double_precision 
  int resMidBitsLength;
  int residualMidBits;
} LossyCompressionElement;
HuffmanTree* createHuffmanTree(int stateNum);
HuffmanTree* createDefaultHuffmanTree();

//_________________________________________ 
extern int  errorBoundMode; 
extern sz_params *confparams_cpr;
//__________________________________________

void symTransform_8bytes(unsigned char data[8]);
void longToBytes_bigEndian(unsigned char *b, uint64_t num);
void symTransform_4bytes(unsigned char data[4]);
void doubleToBytes(unsigned char *b, double num);
unsigned long zlib_compress5(unsigned char* data, unsigned long dataLength, unsigned char** compressBytes, int level);
unsigned long sz_lossless_compress(int losslessCompressor, int level, unsigned char* data, unsigned long dataLength, unsigned char** compressBytes);
void int32ToBytes_bigEndian(unsigned char *b, uint32_t num);
void floatToBytes(unsigned char *b, float num);
void intToBytes_bigEndian(unsigned char *b, unsigned int num);
void int16ToBytes_bigEndian(unsigned char *b, uint16_t num);
void new_DIA(DynamicIntArray **dia, size_t cap) ;
void new_DBA(DynamicByteArray **dba, size_t cap) ;
void addDBA_Data(DynamicByteArray *dba, unsigned char value);
void convertDBAtoBytes(DynamicByteArray *dba, unsigned char** bytes);
void free_DBA(DynamicByteArray *dba);
void free_DIA(DynamicIntArray *dia);
void sizeToBytes(unsigned char* outBytes, size_t size);
int SZ_ReadConf(const char* sz_cfgFile) ;
int SZ_LoadConf(const char* sz_cfgFile) ;
int SZ_Init(const char *configFilePath);
 int getLeftMovingSteps(size_t k, unsigned char resiBitLength);
size_t convertIntArray2ByteArray_fast_dynamic(unsigned char* timeStepType, unsigned char resiBitLength, size_t nbEle, unsigned char **bytes);
size_t convertIntArray2ByteArray_fast_2b(unsigned char* timeStepType, size_t timeStepTypeLength, unsigned char **result);
//-------------------------------------------
node new_node(HuffmanTree* huffmanTree, size_t freq, unsigned int c, node a, node b);
node new_node2(HuffmanTree *huffmanTree, unsigned int c, unsigned char t);
void qinsert(HuffmanTree *huffmanTree, node n);
node qremove(HuffmanTree* huffmanTree);
void build_code(HuffmanTree *huffmanTree, node n, int len, uint64_t out1, uint64_t out2);

void init_new(HuffmanTree* huffmanTree, int *s, size_t length);

void pad_tree_ushort(HuffmanTree* huffmanTree, unsigned short* L, unsigned short* R, unsigned int* C, unsigned char* t, unsigned int i, node root);
void pad_tree_uint(HuffmanTree* huffmanTree, unsigned int* L, unsigned int* R, unsigned int* C, unsigned char* t, unsigned int i, node root);
void pad_tree_uchar(HuffmanTree* huffmanTree, unsigned char* L, unsigned char* R, unsigned int* C, unsigned char* t, unsigned int i, node root);
void SZ_ReleaseHuffman(HuffmanTree* huffmanTree);
unsigned int convert_HuffTree_to_bytes_anyStates(HuffmanTree* huffmanTree, int nodeCount, unsigned char** out) ;
 
//--------------------------- 
void encode(HuffmanTree *huffmanTree, int *s, size_t length, unsigned char *out, size_t *outSize);

void encode_withTree(HuffmanTree* huffmanTree, int *s, size_t length, unsigned char **out, size_t *outSize);
void new_TightDataPointStorageF(TightDataPointStorageF **mythis,
    size_t dataSeriesLength, size_t exactDataNum, 
    int* type, unsigned char* exactMidBytes, size_t exactMidBytes_size,
    unsigned char* leadNumIntArray,  //leadNumIntArray contains readable numbers....
    unsigned char* resiMidBits, size_t resiMidBits_size,
    unsigned char resiBitLength, 
    double realPrecision, float medianValue, char reqLength, unsigned int intervals, 
    unsigned char* pwrErrBoundBytes, size_t pwrErrBoundBytes_size, unsigned char radExpo);

void new_TightDataPointStorageI(TightDataPointStorageI **mythis,
    size_t dataSeriesLength, size_t exactDataNum, int byteSize, 
    int* type, unsigned char* exactDataBytes, size_t exactDataBytes_size,
    double realPrecision, long minValue, int intervals, int dataType); 
    
void new_TightDataPointStorageD(TightDataPointStorageD **mythis, 
    size_t dataSeriesLength, size_t exactDataNum, 
    int* type, unsigned char* exactMidBytes, size_t exactMidBytes_size,
    unsigned char* leadNumIntArray,  //leadNumIntArray contains readable numbers....
    unsigned char* resiMidBits, size_t resiMidBits_size,
    unsigned char resiBitLength, 
    double realPrecision, double medianValue, char reqLength, unsigned int intervals,
    unsigned char* pwrErrBoundBytes, size_t pwrErrBoundBytes_size, unsigned char radExpo);


void addDIA_Data(DynamicIntArray *dia, int value);

void addExactData(DynamicByteArray *exactMidByteArray, DynamicIntArray *exactLeadNumArray, 
    DynamicIntArray *resiBitArray, LossyCompressionElement *lce);
int compIdenticalLeadingBytesCount_float(unsigned char* preBytes, unsigned char* curBytes);
void updateLossyCompElement_Float(unsigned char* curBytes, unsigned char* preBytes, 
    int reqBytesLength, int resiBitsLength,  LossyCompressionElement *lce);

int compIdenticalLeadingBytesCount_double(unsigned char* preBytes, unsigned char* curBytes);

void updateLossyCompElement_Double(unsigned char* curBytes, unsigned char* preBytes, 
    int reqBytesLength, int resiBitsLength,  LossyCompressionElement *lce);

void compressSingleFloatValue(FloatValueCompressElement *vce, float tgtValue, float precision, float medianValue, 
    int reqLength, int reqBytesLength, int resiBitsLength);

void compressInt8Value(int8_t tgtValue, int8_t minValue, int byteSize, unsigned char* bytes);
    
void compressSingleDoubleValue(DoubleValueCompressElement *vce, double tgtValue, double precision, double medianValue, 
int reqLength, int reqBytesLength, int resiBitsLength);

size_t computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
void compressInt16Value(int16_t tgtValue, int16_t minValue, int byteSize, unsigned char* bytes);
void compressInt32Value(int32_t tgtValue, int32_t minValue, int byteSize, unsigned char* bytes);

float computeRangeSize_float(float* oriData, size_t size, float* valueRangeSize, float* medianValue);

double computeRangeSize_double(double* oriData, size_t size, double* valueRangeSize, double* medianValue);
long computeRangeSize_int(void* oriData, int dataType, size_t size, int64_t* valueRangeSize);
  
double getRealPrecision_float(float valueRangeSize, int errBoundMode, double absErrBound, double relBoundRatio, int *status);
double getRealPrecision_int(long valueRangeSize, int errBoundMode, double absErrBound, double relBoundRatio, int *status);
double getRealPrecision_double(double valueRangeSize, int errBoundMode, double absErrBound, double relBoundRatio, int *status);

unsigned int roundUpToPowerOf2(unsigned int base);

void updateQuantizationInfo(int quant_intervals);
short getPrecisionReqLength_double(double precision);

void computeReqLength_float(double realPrecision, short radExpo, int* reqLength, float* medianValue);
void computeReqLength_double(double realPrecision, short radExpo, int* reqLength, double* medianValue);

unsigned int optimize_intervals_float_1D_opt(float *oriData, size_t dataLength, double realPrecision);
unsigned int optimize_intervals_int8_1D(int8_t *oriData, size_t dataLength, double realPrecision);
unsigned int optimize_intervals_int16_1D(int16_t *oriData, size_t dataLength, double realPrecision);
unsigned int optimize_intervals_int32_1D(int32_t *oriData, size_t dataLength, double realPrecision);
unsigned int optimize_intervals_double_1D_opt(double *oriData, size_t dataLength, double realPrecision);

short getExponent_float(float value);

void listAdd_double(double last3CmprsData[3], double value);

short getExponent_double(double value);

TightDataPointStorageF* SZ_compress_float_1D_MDQ(float *oriData, 
size_t dataLength, double realPrecision, float valueRangeSize, float medianValue_f);

TightDataPointStorageD* SZ_compress_double_1D_MDQ(double *oriData, 
size_t dataLength, double realPrecision, double valueRangeSize, double medianValue_d);

void memcpyDBA_Data(DynamicByteArray *dba, unsigned char* data, size_t length);

int computeByteSizePerIntValue(long valueRangeSize);
TightDataPointStorageI* SZ_compress_int8_1D_MDQ(int8_t *oriData, size_t dataLength, double realPrecision, int64_t valueRangeSize, int64_t minValue);
TightDataPointStorageI* SZ_compress_int16_1D_MDQ(int16_t *oriData, size_t dataLength, double realPrecision, int64_t valueRangeSize, int64_t minValue);
TightDataPointStorageI* SZ_compress_int32_1D_MDQ(int32_t *oriData, size_t dataLength, double realPrecision, int64_t valueRangeSize, int64_t minValue);


void convertSZParamsToBytes(sz_params* params, unsigned char* result);

void convertTDPStoBytes_int(TightDataPointStorageI* tdps, unsigned char* bytes, unsigned char sameByte);
void convertTDPStoBytes_float(TightDataPointStorageF* tdps, unsigned char* bytes, unsigned char* dsLengthBytes, unsigned char sameByte);
void convertTDPStoBytes_double(TightDataPointStorageD* tdps, unsigned char* bytes, unsigned char* dsLengthBytes, unsigned char sameByte);

int convertDataTypeSizeCode(int dataTypeSizeCode);
int convertDataTypeSize(int dataTypeSize);

void convertTDPStoFlatBytes_float(TightDataPointStorageF *tdps, unsigned char** bytes, size_t *size);


//Convert TightDataPointStorage to bytes...
void convertTDPStoFlatBytes_int(TightDataPointStorageI *tdps, unsigned char** bytes, size_t *size);
void convertTDPStoFlatBytes_double(TightDataPointStorageD *tdps, unsigned char** bytes, size_t *size);

char SZ_compress_args_float_NoCkRngeNoGzip_1D(unsigned char** newByteData, float *oriData, 
size_t dataLength, double realPrecision, size_t *outSize, float valueRangeSize, float medianValue_f);

char SZ_compress_args_double_NoCkRngeNoGzip_1D(int cmprType, unsigned char** newByteData, double *oriData, 
size_t dataLength, double realPrecision, size_t *outSize, double valueRangeSize, double medianValue_d);

void SZ_compress_args_int8_NoCkRngeNoGzip_1D(unsigned char** newByteData, int8_t *oriData, 
size_t dataLength, double realPrecision, size_t *outSize, int64_t valueRangeSize, int8_t minValue);

void SZ_compress_args_int16_NoCkRngeNoGzip_1D(unsigned char** newByteData, int16_t *oriData, 
size_t dataLength, double realPrecision, size_t *outSize, int64_t valueRangeSize, int16_t minValue);

void SZ_compress_args_int32_NoCkRngeNoGzip_1D(unsigned char** newByteData, int32_t *oriData, 
size_t dataLength, double realPrecision, size_t *outSize, int64_t valueRangeSize, int32_t minValue);

int SZ_compress_args_float(unsigned char** newByteData, float *oriData, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio, double pwRelBoundRatio);

int SZ_compress_args_double(int cmprType, int withRegression, unsigned char** newByteData, double *oriData, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio, double pwRelBoundRatio);

int SZ_compress_args_int8(unsigned char** newByteData, int8_t *oriData, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio);

int SZ_compress_args_int16(unsigned char** newByteData, int16_t *oriData, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio);

int SZ_compress_args_int32(unsigned char** newByteData, int32_t *oriData, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio);

unsigned char* SZ_compress_args(int dataType, void *data, size_t *outSize, int errBoundMode, double absErrBound, 
double relBoundRatio, double pwrBoundRatio, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

double computeABSErrBoundFromPSNR(double psnr, double threshold, double value_range);
unsigned char *SZ_compress(int dataType, void *data, size_t *outSize, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

unsigned long zlib_uncompress5(unsigned char* compressBytes, unsigned long cmpSize, unsigned char** oriData, unsigned long targetOriSize);

unsigned long sz_lossless_decompress(int losslessCompressor, unsigned char* compressBytes, unsigned long cmpSize, unsigned char** oriData, unsigned long targetOriSize);
void new_TightDataPointStorageF_Empty(TightDataPointStorageF **mythis);
void new_TightDataPointStorageI_Empty(TightDataPointStorageI **mythis);
void new_TightDataPointStorageD_Empty(TightDataPointStorageD **mythis);

short bytesToInt16_bigEndian(unsigned char* bytes);
float bytesToFloat(unsigned char* bytes);
int bytesToInt(unsigned char* bytes);
int bytesToInt32_bigEndian(unsigned char* bytes);

double bytesToDouble(unsigned char* bytes);
sz_params* convertBytesToSZParams(unsigned char* bytes);

int bytesToInt_bigEndian(unsigned char* bytes);
long bytesToLong_bigEndian(unsigned char* b);
size_t bytesToSize(unsigned char* bytes);
int new_TightDataPointStorageF_fromFlatBytes(TightDataPointStorageF **mythis, unsigned char* flatBytes, size_t flatBytesLength);
int new_TightDataPointStorageI_fromFlatBytes(TightDataPointStorageI **mythis, unsigned char* flatBytes, size_t flatBytesLength);
int new_TightDataPointStorageD_fromFlatBytes(TightDataPointStorageD **mythis, unsigned char* flatBytes, size_t flatBytesLength); 
int computeDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
void convertByteArray2IntArray_fast_2b(size_t stepLength, unsigned char* byteArray, size_t byteArrayLength, unsigned char **intArray);


void free_TightDataPointStorageD(TightDataPointStorageD *tdps);


void unpad_tree_uchar(HuffmanTree* huffmanTree, unsigned char* L, unsigned char* R, unsigned int* C, unsigned char *t, unsigned int i, node root);

void unpad_tree_ushort(HuffmanTree* huffmanTree, unsigned short* L, unsigned short* R, unsigned int* C, unsigned char* t, unsigned int i, node root);

void unpad_tree_uint(HuffmanTree* huffmanTree, unsigned int* L, unsigned int* R, unsigned int* C, unsigned char* t, unsigned int i, node root);


node reconstruct_HuffTree_from_bytes_anyStates(HuffmanTree *huffmanTree, unsigned char* bytes, int nodeCount);

void decode(unsigned char *s, size_t targetLength, node t, int *out);

void decode_withTree(HuffmanTree* huffmanTree, unsigned char *s, size_t targetLength, int *out);

int getMaskRightCode(int m);
int getRightMovingSteps(int kMod8, int resiBitLength);

int getRightMovingCode(int kMod8, int resiBitLength);
int getLeftMovingCode(int kMod8);

void decompressDataSeries_float_1D(float** data, size_t dataSeriesLength, TightDataPointStorageF* tdps);

void decompressDataSeries_double_1D(double** data, size_t dataSeriesLength, double* hist_data, TightDataPointStorageD* tdps);
int computeRightShiftBits(int exactByteSize, int dataType);
void decompressDataSeries_int8_1D(int8_t** data, size_t dataSeriesLength, TightDataPointStorageI* tdps) ;
void decompressDataSeries_int16_1D(int16_t** data, size_t dataSeriesLength, TightDataPointStorageI* tdps) ;
void decompressDataSeries_int32_1D(int32_t** data, size_t dataSeriesLength, TightDataPointStorageI* tdps);
 


char numberOfLeadingZeros_Int(int i);
char numberOfLeadingZeros_Long(long i) ;
int computeBitNumRequired(size_t dataLength);

void getSnapshotData_float_1D(float** data, size_t dataSeriesLength, TightDataPointStorageF* tdps, int errBoundMode);
void getSnapshotData_int8_1D(int8_t** data, size_t dataSeriesLength, TightDataPointStorageI* tdps, int errBoundMode);
void getSnapshotData_int16_1D(int16_t** data, size_t dataSeriesLength, TightDataPointStorageI* tdps, int errBoundMode);
void getSnapshotData_double_1D(double** data, size_t dataSeriesLength, TightDataPointStorageD* tdps, int errBoundMode, int compressionType, double* hist_data);

 void free_TightDataPointStorageF2(TightDataPointStorageF *tdps);
void free_TightDataPointStorageD2(TightDataPointStorageD *tdps);

int SZ_decompress_args_float(float** newData, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, unsigned char* cmpBytes, size_t cmpSize);
int SZ_decompress_args_int8(int8_t** newData, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, unsigned char* cmpBytes, size_t cmpSize);
int SZ_decompress_args_int16(int16_t** newData, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, unsigned char* cmpBytes, size_t cmpSize);
void getSnapshotData_int32_1D(int32_t** data, size_t dataSeriesLength, TightDataPointStorageI* tdps, int errBoundMode);
int SZ_decompress_args_int32(int32_t** newData, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, unsigned char* cmpBytes, size_t cmpSize);


int SZ_decompress_args_double(double** newData, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, unsigned char* cmpBytes, 
size_t cmpSize, int compressionType, double* hist_data);

void *SZ_decompress(int dataType, unsigned char *bytes, size_t byteLength, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

float rms4float(float *v, int n);
float rms4double(double *v, int n);
float rms4int8(int8_t *v, int n);
float rms4int16(int16_t *v, int n);
float rms4int32(int32_t *v, int n);

#ifdef __cplusplus
}
#endif

#endif
