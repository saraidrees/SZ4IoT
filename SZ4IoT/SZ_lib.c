
#include "SZ_lib.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdint.h>
#include "zlib.h"




int dataEndianType; //*endian type of the data read from disk
int sysEndianType; //*sysEndianType is actually set automatically.

int sz_with_regression;
int MetaDataByteLength=0;
sz_params *confparams_dec;
sz_exedata *exe_params;



void symTransform_8bytes(unsigned char data[8])
{
  unsigned char tmp = data[0];
  data[0] = data[7];
  data[7] = tmp;

  tmp = data[1];
  data[1] = data[6];
  data[6] = tmp;
  
  tmp = data[2];
  data[2] = data[5];
  data[5] = tmp;
  
  tmp = data[3];
  data[3] = data[4];
  data[4] = tmp;
}

void longToBytes_bigEndian(unsigned char *b, uint64_t num)
{
    //printf("%llu \n",num);
//    exit(0);

  b[0] = (unsigned char)(num>>56);
  b[1] = (unsigned char)(num>>48);
  b[2] = (unsigned char)(num>>40);
  b[3] = (unsigned char)(num>>32);
  b[4] = (unsigned char)(num>>24);
  b[5] = (unsigned char)(num>>16);
  b[6] = (unsigned char)(num>>8);
  b[7] = (unsigned char)(num);
//  if(dataEndianType==LITTLE_ENDIAN_DATA)
//    symTransform_8bytes(*b);
}

 void symTransform_4bytes(unsigned char data[4])
{
  unsigned char tmp = data[0];
  data[0] = data[3];
  data[3] = tmp;

  tmp = data[1];
  data[1] = data[2];
  data[2] = tmp;
}


void doubleToBytes(unsigned char *b, double num)
{
  ldouble buf;
  buf.value = num;
  memcpy(b, buf.byte, 8);
  if(sysEndianType==LITTLE_ENDIAN_SYSTEM)
    symTransform_8bytes(b);
}

#define SZ_ZLIB_BUFFER_SIZE 1024


unsigned long zlib_compress5(unsigned char* data, unsigned long dataLength, unsigned char** compressBytes, int level)
{
  int ret, flush;
  unsigned have;
  z_stream strm;
  unsigned char* in = data;

  /* allocate deflate state */
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  ret = deflateInit(&strm, level);
  //int windowBits = 15;
    //ret = deflateInit2(&strm, level, Z_DEFLATED, windowBits, DEF_MEM_LEVEL, Z_DEFAULT_STRATEGY);//Z_FIXED); //Z_DEFAULT_STRATEGY

  if (ret != Z_OK)
    return ret;

  size_t p_size = 0, av_in = 0;
    uLong estCmpLen = deflateBound(&strm, dataLength);
    *compressBytes = (unsigned char*)malloc(sizeof(unsigned char)*estCmpLen); 
  unsigned char* out = *compressBytes; 

  /* compress until end of file */
  do {    
    p_size += SZ_ZLIB_BUFFER_SIZE;
    if(p_size>=dataLength)
    {
      av_in = dataLength - (p_size - SZ_ZLIB_BUFFER_SIZE);
      flush = Z_FINISH;
    }
    else
    {
      av_in = SZ_ZLIB_BUFFER_SIZE;
      flush = Z_NO_FLUSH;
    }
    strm.avail_in = av_in;
    strm.next_in = in;

    /* run deflate() on input until output buffer not full, finish
       compression if all of source has been read in */
    do {
      strm.avail_out = SZ_ZLIB_BUFFER_SIZE;
      strm.next_out = out;
      ret = deflate(&strm, flush);    /* no bad return value */

      have = SZ_ZLIB_BUFFER_SIZE - strm.avail_out;
      out += have;
    } while (strm.avail_out == 0);

    in+=av_in;

    /* done when last data in file processed */
  } while (flush != Z_FINISH);

  /* clean up and return */
  (void)deflateEnd(&strm);  
  
  return strm.total_out;  
}

unsigned long sz_lossless_compress(int losslessCompressor, int level, unsigned char* data, unsigned long dataLength, unsigned char** compressBytes)
{
  unsigned long outSize = 0; 
  size_t estimatedCompressedSize = 0;
  switch(losslessCompressor)
  {
  case GZIP_COMPRESSOR:
 { 
    unsigned long ms1 = micros();
    outSize = zlib_compress5(data, dataLength, compressBytes, level);
    unsigned long ms2 = micros();
   // printf("Time for lossless: %d\n",ms2-ms1);
 }  
    break;

  default:
; //   Serial.println("Error: Unrecognized lossless compressor in sz_lossless_compress()\n");
  }
  return outSize;
}


void int32ToBytes_bigEndian(unsigned char *b, uint32_t num)
{
  b[0] = (unsigned char)(num >> 24);  
  b[1] = (unsigned char)(num >> 16);  
  b[2] = (unsigned char)(num >> 8); 
  b[3] = (unsigned char)(num);    
}

 void floatToBytes(unsigned char *b, float num)
{
  lfloat buf;
  buf.value = num;
  memcpy(b, buf.byte, 4);
  if(sysEndianType==LITTLE_ENDIAN_SYSTEM)
    symTransform_4bytes(b);   
}


void free_TightDataPointStorageF(TightDataPointStorageF *tdps)
{
  if(tdps->rtypeArray!=NULL)
    free(tdps->rtypeArray);
  if(tdps->typeArray!=NULL)
    free(tdps->typeArray);
  if(tdps->leadNumArray!=NULL)
    free(tdps->leadNumArray);
  if(tdps->exactMidBytes!=NULL)
    free(tdps->exactMidBytes);
  if(tdps->residualMidBits!=NULL)
    free(tdps->residualMidBits);
  if(tdps->pwrErrBoundBytes!=NULL)
    free(tdps->pwrErrBoundBytes);
  free(tdps);
}


void free_TightDataPointStorageI(TightDataPointStorageI *tdps)
{
  if(tdps->typeArray!=NULL)
    free(tdps->typeArray);
  if(tdps->exactDataBytes!=NULL)
    free(tdps->exactDataBytes);
  free(tdps);
}

void free_TightDataPointStorageI2(TightDataPointStorageI *tdps)
{
  free(tdps);
}

void free_TightDataPointStorageD(TightDataPointStorageD *tdps)
{
  if(tdps->rtypeArray!=NULL)
    free(tdps->rtypeArray);
  if(tdps->typeArray!=NULL)
    free(tdps->typeArray);
  if(tdps->leadNumArray!=NULL)
    free(tdps->leadNumArray);
  if(tdps->exactMidBytes!=NULL)
    free(tdps->exactMidBytes);
  if(tdps->residualMidBits!=NULL)
    free(tdps->residualMidBits);
  if(tdps->pwrErrBoundBytes!=NULL)  
    free(tdps->pwrErrBoundBytes);
  free(tdps);
}

 void intToBytes_bigEndian(unsigned char *b, unsigned int num)
{
  b[0] = (unsigned char)(num >> 24);  
  b[1] = (unsigned char)(num >> 16);  
  b[2] = (unsigned char)(num >> 8); 
  b[3] = (unsigned char)(num);  
  
  //note: num >> xxx already considered endian_type...
//if(dataEndianType==LITTLE_ENDIAN_DATA)
//    symTransform_4bytes(*b); //change to BIG_ENDIAN_DATA
}

void int16ToBytes_bigEndian(unsigned char *b, uint16_t num)
{
  b[0] = (unsigned char)(num >> 8); 
  b[1] = (unsigned char)(num);
}


void new_DIA(DynamicIntArray **dia, size_t cap) {
    *dia = (DynamicIntArray *)malloc(sizeof(DynamicIntArray));
        (*dia)->size = 0;
        (*dia)->capacity = cap;
        (*dia)->array = (unsigned char*)malloc(sizeof(unsigned char)*cap);
    }



void new_DBA(DynamicByteArray **dba, size_t cap) {
    *dba = (DynamicByteArray *)malloc(sizeof(DynamicByteArray));
        (*dba)->size = 0;
        (*dba)->capacity = cap;
        (*dba)->array = (unsigned char*)malloc(sizeof(unsigned char)*cap);
    }

void addDBA_Data(DynamicByteArray *dba, unsigned char value)
{
  if(dba->size==dba->capacity)
  {
    dba->capacity = dba->capacity << 1;
    dba->array = (unsigned char *)realloc(dba->array, dba->capacity*sizeof(unsigned char));
  }
  dba->array[dba->size] = value;
  dba->size ++;
}


void convertDBAtoBytes(DynamicByteArray *dba, unsigned char** bytes)
{ 
  size_t size = dba->size;
  if(size>0)
    *bytes = (unsigned char*)malloc(size * sizeof(unsigned char));
  else
    *bytes = NULL;
  memcpy(*bytes, dba->array, size*sizeof(unsigned char)); 
}

void free_DBA(DynamicByteArray *dba)
{ 
  free(dba->array);
  free(dba);
}

void free_DIA(DynamicIntArray *dia)
{
  free(dia->array);
  free(dia);
}

void sizeToBytes(unsigned char* outBytes, size_t size)
{ 
  if(exe_params->SZ_SIZE_TYPE==4)
    intToBytes_bigEndian(outBytes, size);//4
  else
    longToBytes_bigEndian(outBytes, size);//8
}


int SZ_ReadConf(const char* sz_cfgFile) { 
    // Check access to SZ configuration file and load dictionary
    //record the setting in confparams_cpr
    confparams_cpr = (sz_params*)malloc(sizeof(sz_params));    
    exe_params = (sz_exedata*)malloc(sizeof(sz_exedata));
    
    int x = 1;
    char sol_name[256];
    char *modeBuf;
    char *errBoundMode;
    char *endianTypeString;
//    dictionary *ini;
    char *par;

  char *y = (char*)&x;
  
  if(*y==1)
    sysEndianType = LITTLE_ENDIAN_SYSTEM;
  else //=0
    sysEndianType = BIG_ENDIAN_SYSTEM;
    
    if(sz_cfgFile == NULL)
    {
    dataEndianType = LITTLE_ENDIAN_DATA;
    confparams_cpr->sol_ID = SZ;
    //RAPH
    confparams_cpr->max_quant_intervals = 256;   //65536
    confparams_cpr->maxRangeRadius = confparams_cpr->max_quant_intervals/2;
        
    exe_params->intvCapacity = confparams_cpr->maxRangeRadius*2;
    exe_params->intvRadius = confparams_cpr->maxRangeRadius;
    
    confparams_cpr->quantization_intervals = 0;
    exe_params->optQuantMode = 1;
    confparams_cpr->predThreshold = 0.99;
    confparams_cpr->sampleDistance = 100;
    
    confparams_cpr->szMode = SZ_BEST_COMPRESSION;
//    confparams_cpr->losslessCompressor = ZSTD_COMPRESSOR; //other option: GZIP_COMPRESSOR;
    confparams_cpr->losslessCompressor = GZIP_COMPRESSOR;
    if(confparams_cpr->losslessCompressor==ZSTD_COMPRESSOR)
      confparams_cpr->gzipMode = 3; //fast mode
    else
      confparams_cpr->gzipMode = 1; //high speed mode
    
    confparams_cpr->errorBoundMode = PSNR;
    confparams_cpr->psnr = 90;
    confparams_cpr->absErrBound = 1E-4;
    confparams_cpr->relBoundRatio = 1E-4;
     
    confparams_cpr->pw_relBoundRatio = 1E-3;
    confparams_cpr->segment_size = 36;
    
    confparams_cpr->pwr_type = SZ_PWR_MIN_TYPE;
    
    confparams_cpr->snapshotCmprStep = 5;
    
    sz_with_regression = SZ_WITH_LINEAR_REGRESSION;
  
    confparams_cpr->randomAccess = 0; //0: no random access , 1: support random access
  
    return SZ_SCES;
  }
    
  
    return SZ_SCES;
}



int SZ_LoadConf(const char* sz_cfgFile) { 
    int res = SZ_ReadConf(sz_cfgFile);
    if (res != SZ_SCES)
    {
       // printf("[SZ] ERROR: Impossible to read configuration.\n");
        return SZ_NSCS;
    }
    return SZ_SCES;
}


int SZ_Init(const char *configFilePath)
{ 
  int loadFileResult = SZ_LoadConf(configFilePath);
  if(loadFileResult==SZ_NSCS)
    return SZ_NSCS;
  
  exe_params->SZ_SIZE_TYPE = sizeof(size_t);
  

  return SZ_SCES;
}


 int getLeftMovingSteps(size_t k, unsigned char resiBitLength)
{ 
  return 8 - k%8 - resiBitLength;
}



size_t convertIntArray2ByteArray_fast_dynamic(unsigned char* timeStepType, unsigned char resiBitLength, size_t nbEle, unsigned char **bytes)
{ 
  size_t i = 0, j = 0, k = 0; 
  int value;
  DynamicByteArray* dba;
  new_DBA(&dba, 1024);
  int tmp = 0, leftMovSteps = 0;
  for(j = 0;j<nbEle;j++)
  {
    if(resiBitLength==0)
      continue;
    value = timeStepType[i];
    leftMovSteps = getLeftMovingSteps(k, resiBitLength);
    if(leftMovSteps < 0)
    {
      tmp = tmp | (value >> (-leftMovSteps));
      addDBA_Data(dba, (unsigned char)tmp);
      tmp = 0 | (value << (8+leftMovSteps));
    }
    else if(leftMovSteps > 0)
    {
      tmp = tmp | (value << leftMovSteps);
    }
    else //==0
    {
      tmp = tmp | value;
      addDBA_Data(dba, (unsigned char)tmp);
      tmp = 0;
    }
    i++;
    k += resiBitLength;
  }
  if(leftMovSteps != 0)
    addDBA_Data(dba, (unsigned char)tmp);
  convertDBAtoBytes(dba, bytes);
  size_t size = dba->size;
  free_DBA(dba);
  return size;
}


size_t convertIntArray2ByteArray_fast_2b(unsigned char* timeStepType, size_t timeStepTypeLength, unsigned char **result)
{
  size_t i, j, byteLength = 0;
  if(timeStepTypeLength%4==0)
    byteLength = timeStepTypeLength*2/8;
  else
    byteLength = timeStepTypeLength*2/8+1;
  if(byteLength>0)
    *result = (unsigned char*)malloc(byteLength*sizeof(unsigned char));
  else
    *result = NULL;
  size_t n = 0;
  for(i = 0;i<byteLength;i++)
  {
    int tmp = 0;
    for(j = 0;j<4&&n<timeStepTypeLength;j++)
    {
      int type = timeStepType[n];
      switch(type)
      {
      case 0: 
        
        break;
      case 1:
        tmp = (tmp | (1 << (6-j*2)));
        break;
      case 2:
        tmp = (tmp | (2 << (6-j*2)));
        break;
      case 3:
        tmp = (tmp | (3 << (6-j*2)));
        break;
      default:
       // printf("Error: wrong timestep type...: type[%zu]=%d\n", n, type);
        exit(0);
      }
      n++;
    }
    (*result)[i] = (unsigned char)tmp;
  }
  return byteLength;
}


void pad_tree_ushort(HuffmanTree* huffmanTree, unsigned short* L, unsigned short* R, unsigned int* C, unsigned char* t, unsigned int i, node root)
{
  C[i] = root->c;
  t[i] = root->t;
  node lroot = root->left;
  if(lroot!=0)
  {
    huffmanTree->n_inode++;
    L[i] = huffmanTree->n_inode;
    pad_tree_ushort(huffmanTree,L,R,C,t,huffmanTree->n_inode, lroot);
  }
  node rroot = root->right;
  if(rroot!=0)
  {
    huffmanTree->n_inode++;
    R[i] = huffmanTree->n_inode;
    pad_tree_ushort(huffmanTree,L,R,C,t,huffmanTree->n_inode, rroot);
  } 
}

void pad_tree_uint(HuffmanTree* huffmanTree, unsigned int* L, unsigned int* R, unsigned int* C, unsigned char* t, unsigned int i, node root)
{
  C[i] = root->c;
  t[i] = root->t;
  node lroot = root->left;
  if(lroot!=0)
  {
    huffmanTree->n_inode++;
    L[i] = huffmanTree->n_inode;
    pad_tree_uint(huffmanTree,L,R,C,t,huffmanTree->n_inode, lroot);
  }
  node rroot = root->right;
  if(rroot!=0)
  {
    huffmanTree->n_inode++;
    R[i] = huffmanTree->n_inode;
    pad_tree_uint(huffmanTree,L,R,C,t,huffmanTree->n_inode, rroot);
  }
}

void pad_tree_uchar(HuffmanTree* huffmanTree, unsigned char* L, unsigned char* R, unsigned int* C, unsigned char* t, unsigned int i, node root)
{
  C[i] = root->c;
  t[i] = root->t;
  node lroot = root->left;
  if(lroot!=0)
  {
    huffmanTree->n_inode++;
    L[i] = huffmanTree->n_inode;
    pad_tree_uchar(huffmanTree, L,R,C,t, huffmanTree->n_inode, lroot);
  }
  node rroot = root->right;
  if(rroot!=0)
  {
    huffmanTree->n_inode++;
    R[i] = huffmanTree->n_inode;
    pad_tree_uchar(huffmanTree, L,R,C,t, huffmanTree->n_inode, rroot);
  }
}
 

  
void SZ_ReleaseHuffman(HuffmanTree* huffmanTree)
{
  size_t i;
  free(huffmanTree->pool);
  huffmanTree->pool = NULL;
  free(huffmanTree->qqq);
  huffmanTree->qqq = NULL;
  for(i=0;i<huffmanTree->stateNum;i++)
  {
    if(huffmanTree->code[i]!=NULL)
      free(huffmanTree->code[i]);
  }
  free(huffmanTree->code);
  huffmanTree->code = NULL;
  free(huffmanTree->cout);
  huffmanTree->cout = NULL; 
  free(huffmanTree);
  huffmanTree = NULL;
}


unsigned int convert_HuffTree_to_bytes_anyStates(HuffmanTree* huffmanTree, int nodeCount, unsigned char** out) 
{
  //printf("nodeCount=%d\n", nodeCount);
  if(nodeCount<=256)
  {
    unsigned char* L = (unsigned char*)malloc(nodeCount*sizeof(unsigned char));
    memset(L, 0, nodeCount*sizeof(unsigned char));
    unsigned char* R = (unsigned char*)malloc(nodeCount*sizeof(unsigned char));
    memset(R, 0, nodeCount*sizeof(unsigned char));
    unsigned int* C = (unsigned int*)malloc(nodeCount*sizeof(unsigned int));
    memset(C, 0, nodeCount*sizeof(unsigned int));
    unsigned char* t = (unsigned char*)malloc(nodeCount*sizeof(unsigned char));
    memset(t, 0, nodeCount*sizeof(unsigned char));

    pad_tree_uchar(huffmanTree,L,R,C,t,0,huffmanTree->qq[1]);

    unsigned int totalSize = 1+3*nodeCount*sizeof(unsigned char)+nodeCount*sizeof(unsigned int);  
    *out = (unsigned char*)malloc(totalSize*sizeof(unsigned char));
    (*out)[0] = (unsigned char)sysEndianType;
    memcpy(*out+1, L, nodeCount*sizeof(unsigned char));
    memcpy((*out)+1+nodeCount*sizeof(unsigned char),R,nodeCount*sizeof(unsigned char));
    memcpy((*out)+1+2*nodeCount*sizeof(unsigned char),C,nodeCount*sizeof(unsigned int));
    memcpy((*out)+1+2*nodeCount*sizeof(unsigned char)+nodeCount*sizeof(unsigned int), t, nodeCount*sizeof(unsigned char));
    free(L);
    free(R);
    free(C);
    free(t);
    return totalSize;

  }
  else if(nodeCount<=65536)
  {
    unsigned short* L = (unsigned short*)malloc(nodeCount*sizeof(unsigned short));
    memset(L, 0, nodeCount*sizeof(unsigned short));
    unsigned short* R = (unsigned short*)malloc(nodeCount*sizeof(unsigned short));
    memset(R, 0, nodeCount*sizeof(unsigned short));
    unsigned int* C = (unsigned int*)malloc(nodeCount*sizeof(unsigned int));  
    memset(C, 0, nodeCount*sizeof(unsigned int));   
    unsigned char* t = (unsigned char*)malloc(nodeCount*sizeof(unsigned char));
    memset(t, 0, nodeCount*sizeof(unsigned char));    
    pad_tree_ushort(huffmanTree,L,R,C,t,0,huffmanTree->qq[1]);
    unsigned int totalSize = 1+2*nodeCount*sizeof(unsigned short)+nodeCount*sizeof(unsigned char) + nodeCount*sizeof(unsigned int);
    *out = (unsigned char*)malloc(totalSize);
    (*out)[0] = (unsigned char)sysEndianType;   
    memcpy(*out+1, L, nodeCount*sizeof(unsigned short));
    memcpy((*out)+1+nodeCount*sizeof(unsigned short),R,nodeCount*sizeof(unsigned short));
    memcpy((*out)+1+2*nodeCount*sizeof(unsigned short),C,nodeCount*sizeof(unsigned int));
    memcpy((*out)+1+2*nodeCount*sizeof(unsigned short)+nodeCount*sizeof(unsigned int),t,nodeCount*sizeof(unsigned char));
    free(L);
    free(R);
    free(C);
    free(t);    
    return totalSize;
  }
  else //nodeCount>65536
  {
    unsigned int* L = (unsigned int*)malloc(nodeCount*sizeof(unsigned int));
    memset(L, 0, nodeCount*sizeof(unsigned int));
    unsigned int* R = (unsigned int*)malloc(nodeCount*sizeof(unsigned int));
    memset(R, 0, nodeCount*sizeof(unsigned int));
    unsigned int* C = (unsigned int*)malloc(nodeCount*sizeof(unsigned int));  
    memset(C, 0, nodeCount*sizeof(unsigned int));
    unsigned char* t = (unsigned char*)malloc(nodeCount*sizeof(unsigned char));
    memset(t, 0, nodeCount*sizeof(unsigned char));
    pad_tree_uint(huffmanTree, L,R,C,t,0,huffmanTree->qq[1]);
    
    //debug
    //node root = new_node2(0,0);
    //unpad_tree_uint(L,R,C,t,0,root);    
    
    unsigned int totalSize = 1+3*nodeCount*sizeof(unsigned int)+nodeCount*sizeof(unsigned char);
    *out = (unsigned char*)malloc(totalSize);
    (*out)[0] = (unsigned char)sysEndianType;
    memcpy(*out+1, L, nodeCount*sizeof(unsigned int));
    memcpy((*out)+1+nodeCount*sizeof(unsigned int),R,nodeCount*sizeof(unsigned int));
    memcpy((*out)+1+2*nodeCount*sizeof(unsigned int),C,nodeCount*sizeof(unsigned int));
    memcpy((*out)+1+3*nodeCount*sizeof(unsigned int),t,nodeCount*sizeof(unsigned char));
    free(L);
    free(R);
    free(C);
    free(t);
    return totalSize;   
  }
}

node new_node(HuffmanTree* huffmanTree, size_t freq, unsigned int c, node a, node b)
{
  node n = huffmanTree->pool + huffmanTree->n_nodes++;
  if (freq) 
  {
    n->c = c;
    n->freq = freq;
    n->t = 1;
  }
  else {
    n->left = a; 
    n->right = b;
    n->freq = a->freq + b->freq;
    n->t = 0;
    //n->c = 0;
  }
  return n;
}
 
node qremove(HuffmanTree* huffmanTree)
{ 
  int i, l;
  node n = huffmanTree->qq[i = 1];
 
  if (huffmanTree->qend < 2) return 0;
  huffmanTree->qend --;
  while ((l = (i<<1)) < huffmanTree->qend)  //l=(i*2)
  {
    if (l + 1 < huffmanTree->qend && huffmanTree->qq[l + 1]->freq < huffmanTree->qq[l]->freq) l++;
    huffmanTree->qq[i] = huffmanTree->qq[l], i = l;
  }
  huffmanTree->qq[i] = huffmanTree->qq[huffmanTree->qend];
  return n;
}
void build_code(HuffmanTree *huffmanTree, node n, int len, uint64_t out1, uint64_t out2)
{
  if (n->t) {
      //printf("hufman build n->c %d out1 %u out2 %u \n",n->c,out1,out2);
      //exit(0);
    huffmanTree->code[n->c] = (uint64_t*)malloc(2*sizeof(uint64_t));
    if(len<=64)
    {
      (huffmanTree->code[n->c])[0] = out1 << (64 - len);
      (huffmanTree->code[n->c])[1] = out2;
    }
    else
    {
      (huffmanTree->code[n->c])[0] = out1;
      (huffmanTree->code[n->c])[1] = out2 << (128 - len);
    }
    huffmanTree->cout[n->c] = (unsigned char)len;
    return;
  }
  int index = len >> 6; //=len/64
  if(index == 0)
  {
    out1 = out1 << 1;
    out1 = out1 | 0;
    build_code(huffmanTree, n->left, len + 1, out1, 0);
    out1 = out1 | 1;
    build_code(huffmanTree, n->right, len + 1, out1, 0);    
  }
  else
  {
    if(len%64!=0)
      out2 = out2 << 1;
    out2 = out2 | 0;
    build_code(huffmanTree, n->left, len + 1, out1, out2);
    out2 = out2 | 1;
    build_code(huffmanTree, n->right, len + 1, out1, out2); 
  }
}


void qinsert(HuffmanTree *huffmanTree, node n)
{
  int j, i = huffmanTree->qend++;
  while ((j = (i>>1)))  //j=i/2
  {
    if (huffmanTree->qq[j]->freq <= n->freq) break;
    huffmanTree->qq[i] = huffmanTree->qq[j], i = j;
  }
  huffmanTree->qq[i] = n;
}

void init_new(HuffmanTree* huffmanTree, int *s, size_t length)
{
  size_t i, index;
  size_t *freq = (size_t *)malloc(huffmanTree->allNodes*sizeof(size_t));
  memset(freq, 0, huffmanTree->allNodes*sizeof(size_t));
  for(i = 0;i < length;i++) 
  {
    //index = 0;
    //index = (index | s[i])<<8;
    //index = index | s[i+1];
    index = s[i];
    freq[index]++;
  }
 
  for (i = 0; i < huffmanTree->allNodes; i++)
    if (freq[i]) 
      qinsert(huffmanTree, new_node(huffmanTree, freq[i], i, 0, 0));
 
  while (huffmanTree->qend > 2) 
    qinsert(huffmanTree, new_node(huffmanTree, 0, 0, qremove(huffmanTree), qremove(huffmanTree)));
 
  build_code(huffmanTree, huffmanTree->qq[1], 0, 0, 0);
  free(freq);
}
 
void encode(HuffmanTree *huffmanTree, int *s, size_t length, unsigned char *out, size_t *outSize)
{
  size_t i = 0;
  unsigned char bitSize = 0, byteSize, byteSizep;
  int state;
  unsigned char *p = out;
  int lackBits = 0;
  //long totalBitSize = 0, maxBitSize = 0, bitSize21 = 0, bitSize32 = 0;
  for (i = 0;i<length;i++) 
  {
    //state = 0;
    //state = (state | s[i])<<8;
    //state = state | s[i+1];
    
    state = s[i];
    bitSize = huffmanTree->cout[state]; 
    
    //printf("%d %d : %d %u\n",i, state, bitSize, (code[state])[0] >> (64-cout[state])); 
    //debug: compute the average bitSize and the count that is over 32...   
    /*if(bitSize>=21)
      bitSize21++;
    if(bitSize>=32)
      bitSize32++;
    if(maxBitSize<bitSize)
      maxBitSize = bitSize;
    totalBitSize+=bitSize;*/

    if(lackBits==0)
    {
      byteSize = bitSize%8==0 ? bitSize/8 : bitSize/8+1; //it's equal to the number of bytes involved (for *outSize)
      byteSizep = bitSize/8; //it's used to move the pointer p for next data
      if(byteSize<=8)       
      {
        longToBytes_bigEndian(p, (huffmanTree->code[state])[0]);
        p += byteSizep;
      }
      else //byteSize>8
      {
        longToBytes_bigEndian(p, (huffmanTree->code[state])[0]);
        p += 8;     
        longToBytes_bigEndian(p, (huffmanTree->code[state])[1]);
        p += (byteSizep - 8);   
      }
      *outSize += byteSize;
      lackBits = bitSize%8==0 ? 0 : 8 - bitSize%8;
    }
    else
    {
        //printf("code %lu\n",huffmanTree->code[state]);
      *p = (*p) | (unsigned char)((huffmanTree->code[state])[0] >> (64 - lackBits));      
      if(lackBits < bitSize)
      {
        p++;
        //(*outSize)++;
        uint64_t newCode = (huffmanTree->code[state])[0] << lackBits;
        longToBytes_bigEndian(p, newCode);        

        if(bitSize<=64)
        {
          bitSize -= lackBits;
          byteSize = bitSize%8==0 ? bitSize/8 : bitSize/8+1;
          byteSizep = bitSize/8;
          p += byteSizep;
          (*outSize)+=byteSize;
          lackBits = bitSize%8==0 ? 0 : 8 - bitSize%8;
        }
        else //bitSize > 64
        {
          byteSizep = 7; //must be 7 bytes, because lackBits!=0
          p+=byteSizep;
          (*outSize)+=byteSize;
          
          bitSize -= 64;
          if(lackBits < bitSize)
          {
            *p = (*p) | (unsigned char)((huffmanTree->code[state])[0] >> (64 - lackBits));
            p++;
            //(*outSize)++;           
            newCode = (huffmanTree->code[state])[1] << lackBits;
            longToBytes_bigEndian(p, newCode);
            bitSize -= lackBits;
            byteSize = bitSize%8==0 ? bitSize/8 : bitSize/8+1;
            byteSizep = bitSize/8;
            p += byteSizep;
            (*outSize)+=byteSize;
            lackBits = bitSize%8==0 ? 0 : 8 - bitSize%8;            
          }
          else //lackBits >= bitSize
          {
            *p = (*p) | (unsigned char)((huffmanTree->code[state])[0] >> (64 - bitSize));
            lackBits -= bitSize;
          }   
        }
      }
      else //lackBits >= bitSize
      {
        lackBits -= bitSize;
        if(lackBits==0)
          p++;
      }
    }
  }

 }



void encode_withTree(HuffmanTree* huffmanTree, int *s, size_t length, unsigned char **out, size_t *outSize)
{
  size_t i; 
  int nodeCount = 0;
  unsigned char *treeBytes, buffer[4];
  
  init_new(huffmanTree, s, length);
  for (i = 0; i < huffmanTree->stateNum; i++)
    if (huffmanTree->code[i]) nodeCount++; 
  nodeCount = nodeCount*2-1;
  unsigned int treeByteSize = convert_HuffTree_to_bytes_anyStates(huffmanTree,nodeCount, &treeBytes);
  //printf("treeByteSize=%d\n", treeByteSize);
  *out = (unsigned char*)malloc(length*sizeof(int)+treeByteSize);
  intToBytes_bigEndian(buffer, nodeCount);
  memcpy(*out, buffer, 4);
  intToBytes_bigEndian(buffer, huffmanTree->stateNum/2); //real number of intervals
  memcpy(*out+4, buffer, 4);
  memcpy(*out+8, treeBytes, treeByteSize);
  free(treeBytes);
  size_t enCodeSize = 0;
  encode(huffmanTree, s, length, *out+8+treeByteSize, &enCodeSize);
  *outSize = 8+treeByteSize+enCodeSize;
  
}





HuffmanTree* createHuffmanTree(int stateNum)
{   
  HuffmanTree *huffmanTree = (HuffmanTree*)malloc(sizeof(HuffmanTree));
  memset(huffmanTree, 0, sizeof(HuffmanTree));
  huffmanTree->stateNum = stateNum;
  huffmanTree->allNodes = 2*stateNum;
  
  huffmanTree->pool = (struct node_t*)malloc(huffmanTree->allNodes*2*sizeof(struct node_t));
  huffmanTree->qqq = (node*)malloc(huffmanTree->allNodes*2*sizeof(node));
  huffmanTree->code = (uint64_t**)malloc(huffmanTree->stateNum*sizeof(uint64_t*));
  huffmanTree->cout = (unsigned char *)malloc(huffmanTree->stateNum*sizeof(unsigned char));
  
  memset(huffmanTree->pool, 0, huffmanTree->allNodes*2*sizeof(struct node_t));
  memset(huffmanTree->qqq, 0, huffmanTree->allNodes*2*sizeof(node));
  //printf("size huffman code %d\n",huffmanTree->stateNum*sizeof(uint64_t*));
    memset(huffmanTree->code, 0, huffmanTree->stateNum*sizeof(uint64_t*));
    //for(int i=0;i<huffmanTree->stateNum;i++)
    //    huffmanTree->code[i]=0;


    memset(huffmanTree->cout, 0, huffmanTree->stateNum*sizeof(unsigned char));

  huffmanTree->qq = huffmanTree->qqq - 1;
  huffmanTree->n_nodes = 0;
    huffmanTree->n_inode = 0;
    huffmanTree->qend = 1;  
    
    return huffmanTree;
}

void new_TightDataPointStorageF(TightDataPointStorageF **mythis,
    size_t dataSeriesLength, size_t exactDataNum, 
    int* type, unsigned char* exactMidBytes, size_t exactMidBytes_size,
    unsigned char* leadNumIntArray,  //leadNumIntArray contains readable numbers....
    unsigned char* resiMidBits, size_t resiMidBits_size,
    unsigned char resiBitLength, 
    double realPrecision, float medianValue, char reqLength, unsigned int intervals, 
    unsigned char* pwrErrBoundBytes, size_t pwrErrBoundBytes_size, unsigned char radExpo) {  
  
  *mythis = (TightDataPointStorageF *)malloc(sizeof(TightDataPointStorageF));
  (*mythis)->allSameData = 0;
  (*mythis)->realPrecision = realPrecision;
  (*mythis)->medianValue = medianValue;
  (*mythis)->reqLength = reqLength;

  (*mythis)->dataSeriesLength = dataSeriesLength;
  (*mythis)->exactDataNum = exactDataNum;

  (*mythis)->rtypeArray = NULL;
  (*mythis)->rtypeArray_size = 0;

  int stateNum = 2*intervals;
  HuffmanTree* huffmanTree = createHuffmanTree(stateNum);
  encode_withTree(huffmanTree, type, dataSeriesLength, &(*mythis)->typeArray, &(*mythis)->typeArray_size);
  SZ_ReleaseHuffman(huffmanTree);
    
  (*mythis)->exactMidBytes = exactMidBytes;
  (*mythis)->exactMidBytes_size = exactMidBytes_size;

  (*mythis)->leadNumArray_size = convertIntArray2ByteArray_fast_2b(leadNumIntArray, exactDataNum, &((*mythis)->leadNumArray));

  (*mythis)->residualMidBits_size = convertIntArray2ByteArray_fast_dynamic(resiMidBits, resiBitLength, exactDataNum, &((*mythis)->residualMidBits));
  
  (*mythis)->intervals = intervals;
  
  (*mythis)->isLossless = 0;
  
    
    (*mythis)->pwrErrBoundBytes = NULL;
    
  (*mythis)->radExpo = radExpo;
  
  (*mythis)->pwrErrBoundBytes_size = pwrErrBoundBytes_size;
}




void new_TightDataPointStorageI(TightDataPointStorageI **mythis,
    size_t dataSeriesLength, size_t exactDataNum, int byteSize, 
    int* type, unsigned char* exactDataBytes, size_t exactDataBytes_size,
    double realPrecision, long minValue, int intervals, int dataType) 
{

  *mythis = (TightDataPointStorageI *)malloc(sizeof(TightDataPointStorageI));
  
  //int i = 0;
  (*mythis)->allSameData = 0;
  (*mythis)->realPrecision = realPrecision;
  (*mythis)->minValue = minValue;
  (*mythis)->dataTypeSize = 1;


  (*mythis)->dataSeriesLength = dataSeriesLength;
  (*mythis)->exactDataNum = exactDataNum;
  (*mythis)->exactByteSize = byteSize;


  int stateNum = 2*intervals;
  HuffmanTree* huffmanTree = createHuffmanTree(stateNum);
  encode_withTree(huffmanTree, type, dataSeriesLength, &(*mythis)->typeArray, &(*mythis)->typeArray_size);
  SZ_ReleaseHuffman(huffmanTree);
    
  (*mythis)->exactDataBytes = exactDataBytes;
  (*mythis)->exactDataBytes_size = exactDataBytes_size;
  
  (*mythis)->intervals = intervals;
  
  (*mythis)->isLossless = 0;
}





void new_TightDataPointStorageD(TightDataPointStorageD **mythis, 
    size_t dataSeriesLength, size_t exactDataNum, 
    int* type, unsigned char* exactMidBytes, size_t exactMidBytes_size,
    unsigned char* leadNumIntArray,  //leadNumIntArray contains readable numbers....
    unsigned char* resiMidBits, size_t resiMidBits_size,
    unsigned char resiBitLength, 
    double realPrecision, double medianValue, char reqLength, unsigned int intervals,
    unsigned char* pwrErrBoundBytes, size_t pwrErrBoundBytes_size, unsigned char radExpo) {
  //int i = 0;
  *mythis = (TightDataPointStorageD *)malloc(sizeof(TightDataPointStorageD));
  (*mythis)->allSameData = 0;
  (*mythis)->realPrecision = realPrecision;
  (*mythis)->medianValue = medianValue;
  (*mythis)->reqLength = reqLength;

  (*mythis)->dataSeriesLength = dataSeriesLength;
  (*mythis)->exactDataNum = exactDataNum;

  (*mythis)->rtypeArray = NULL;
  (*mythis)->rtypeArray_size = 0;

  int stateNum = 2*intervals;
  HuffmanTree* huffmanTree = createHuffmanTree(stateNum);
//sara if(confparams_cpr->errorBoundMode == PW_REL && confparams_cpr->accelerate_pw_rel_compression)
 //sara   (*mythis)->max_bits = encode_withTree_MSST19(huffmanTree, type, dataSeriesLength, &(*mythis)->typeArray, &(*mythis)->typeArray_size);
 // else
    encode_withTree(huffmanTree, type, dataSeriesLength, &(*mythis)->typeArray, &(*mythis)->typeArray_size);
  SZ_ReleaseHuffman(huffmanTree);
    
  (*mythis)->exactMidBytes = exactMidBytes;
  (*mythis)->exactMidBytes_size = exactMidBytes_size;

  (*mythis)->leadNumArray_size = convertIntArray2ByteArray_fast_2b(leadNumIntArray, exactDataNum, &((*mythis)->leadNumArray));

  (*mythis)->residualMidBits_size = convertIntArray2ByteArray_fast_dynamic(resiMidBits, resiBitLength, exactDataNum, &((*mythis)->residualMidBits));
  
  (*mythis)->intervals = intervals;
  
  (*mythis)->isLossless = 0;
  
  if(confparams_cpr->errorBoundMode>=PW_REL)
    (*mythis)->pwrErrBoundBytes = pwrErrBoundBytes;
  else
    (*mythis)->pwrErrBoundBytes = NULL;
    
  (*mythis)->radExpo = radExpo;
  
  (*mythis)->pwrErrBoundBytes_size = pwrErrBoundBytes_size;
}


 void listAdd_float(float last3CmprsData[3], float value)
{
  last3CmprsData[2] = last3CmprsData[1];
  last3CmprsData[1] = last3CmprsData[0];
  last3CmprsData[0] = value;
}

 void listAdd_int(int64_t last3CmprsData[3], int64_t value)
{
  last3CmprsData[2] = last3CmprsData[1];
  last3CmprsData[1] = last3CmprsData[0];
  last3CmprsData[0] = value;
}

 void listAdd_double(double last3CmprsData[3], double value)
{
  last3CmprsData[2] = last3CmprsData[1];
  last3CmprsData[1] = last3CmprsData[0];
  last3CmprsData[0] = value;
}

void addDIA_Data(DynamicIntArray *dia, int value)
{
  if(dia->size==dia->capacity)
  {
    dia->capacity = dia->capacity << 1;
    dia->array = (unsigned char *)realloc(dia->array, dia->capacity*sizeof(unsigned char));
  }
  dia->array[dia->size] = (unsigned char)value;
  dia->size ++;
}


void addExactData(DynamicByteArray *exactMidByteArray, DynamicIntArray *exactLeadNumArray, 
    DynamicIntArray *resiBitArray, LossyCompressionElement *lce)
{
  int i;
  int leadByteLength = lce->leadingZeroBytes;
  addDIA_Data(exactLeadNumArray, leadByteLength);
  unsigned char* intMidBytes = lce->integerMidBytes;
  int integerMidBytesLength = lce->integerMidBytes_Length;
  int resMidBitsLength = lce->resMidBitsLength;
  if(intMidBytes!=NULL||resMidBitsLength!=0)
  {
    if(intMidBytes!=NULL)
      for(i = 0;i<integerMidBytesLength;i++)
        addDBA_Data(exactMidByteArray, intMidBytes[i]);
    if(resMidBitsLength!=0)
      addDIA_Data(resiBitArray, lce->residualMidBits);
  }
}


int compIdenticalLeadingBytesCount_float(unsigned char* preBytes, unsigned char* curBytes)
{
  int i, n = 0;
  for(i=0;i<4;i++)
    if(preBytes[i]==curBytes[i])
      n++;
    else
      break;
  if(n>3) n = 3;
  return n;
}

int compIdenticalLeadingBytesCount_double(unsigned char* preBytes, unsigned char* curBytes)
{
  int i, n = 0;
  for(i=0;i<8;i++)
    if(preBytes[i]==curBytes[i])
      n++;
    else
      break;
  if(n>3) n = 3;
  return n;
}


void updateLossyCompElement_Float(unsigned char* curBytes, unsigned char* preBytes, 
    int reqBytesLength, int resiBitsLength,  LossyCompressionElement *lce)
{
  int resiIndex, intMidBytes_Length = 0;
  int leadingNum = compIdenticalLeadingBytesCount_float(preBytes, curBytes); //in fact, float is enough for both single-precision and double-precisiond ata.
  int fromByteIndex = leadingNum;
  int toByteIndex = reqBytesLength; //later on: should use "< toByteIndex" to tarverse....
  if(fromByteIndex < toByteIndex)
  {
    intMidBytes_Length = reqBytesLength - leadingNum;
    memcpy(lce->integerMidBytes, &(curBytes[fromByteIndex]), intMidBytes_Length);
  }
  int resiBits = 0;
  if(resiBitsLength!=0)
  {
    resiIndex = reqBytesLength;
    if(resiIndex < 8)
      resiBits = (curBytes[resiIndex] & 0xFF) >> (8-resiBitsLength);
  }
  lce->leadingZeroBytes = leadingNum;
  lce->integerMidBytes_Length = intMidBytes_Length;
  lce->resMidBitsLength = resiBitsLength;
  lce->residualMidBits = resiBits;
}

void updateLossyCompElement_Double(unsigned char* curBytes, unsigned char* preBytes, 
    int reqBytesLength, int resiBitsLength,  LossyCompressionElement *lce)
{
  int resiIndex, intMidBytes_Length = 0;
  int leadingNum = compIdenticalLeadingBytesCount_double(preBytes, curBytes); //in fact, float is enough for both single-precision and double-precisiond ata.
  int fromByteIndex = leadingNum;
  int toByteIndex = reqBytesLength; //later on: should use "< toByteIndex" to tarverse....
  if(fromByteIndex < toByteIndex)
  {
    intMidBytes_Length = reqBytesLength - leadingNum;
    memcpy(lce->integerMidBytes, &(curBytes[fromByteIndex]), intMidBytes_Length);
  }
  int resiBits = 0;
  if(resiBitsLength!=0)
  {
    resiIndex = reqBytesLength;
    if(resiIndex < 8)
      resiBits = (curBytes[resiIndex] & 0xFF) >> (8-resiBitsLength);
  }
  lce->leadingZeroBytes = leadingNum;
  lce->integerMidBytes_Length = intMidBytes_Length;
  lce->resMidBitsLength = resiBitsLength;
  lce->residualMidBits = resiBits;
}


void compressSingleFloatValue(FloatValueCompressElement *vce, float tgtValue, float precision, float medianValue, 
    int reqLength, int reqBytesLength, int resiBitsLength)
{   
  float normValue = tgtValue - medianValue;

  lfloat lfBuf;
  lfBuf.value = normValue;
      
  int ignBytesLength = 32 - reqLength;
  if(ignBytesLength<0)
    ignBytesLength = 0;
  
  int tmp_int = lfBuf.ivalue;
  intToBytes_bigEndian(vce->curBytes, tmp_int);
    
  lfBuf.ivalue = (lfBuf.ivalue >> ignBytesLength) << ignBytesLength;
  
  //float tmpValue = lfBuf.value;
  
  vce->data = lfBuf.value+medianValue;
  vce->curValue = tmp_int;
  vce->reqBytesLength = reqBytesLength;
  vce->resiBitsLength = resiBitsLength;
}


 void compressInt8Value(int8_t tgtValue, int8_t minValue, int byteSize, unsigned char* bytes)
{
  uint8_t data = tgtValue - minValue;
  memcpy(bytes, &data, byteSize); //byteSize==1
}

 void compressInt16Value(int16_t tgtValue, int16_t minValue, int byteSize, unsigned char* bytes)
{
  uint16_t data = tgtValue - minValue;
  unsigned char tmpBytes[2];
  int16ToBytes_bigEndian(tmpBytes, data);
  memcpy(bytes, tmpBytes + 2 - byteSize, byteSize);
}

 void compressInt32Value(int32_t tgtValue, int32_t minValue, int byteSize, unsigned char* bytes)
{
  uint32_t data = tgtValue - minValue;
  unsigned char tmpBytes[4];
  int32ToBytes_bigEndian(tmpBytes, data);
  memcpy(bytes, tmpBytes + 4 - byteSize, byteSize);
}


void compressSingleDoubleValue(DoubleValueCompressElement *vce, double tgtValue, double precision, double medianValue, 
    int reqLength, int reqBytesLength, int resiBitsLength)
{   
  double normValue = tgtValue - medianValue;

  ldouble lfBuf;
  lfBuf.value = normValue;
      
  int ignBytesLength = 64 - reqLength;
  if(ignBytesLength<0)
    ignBytesLength = 0;

  long long tmp_long = lfBuf.lvalue;
  longToBytes_bigEndian(vce->curBytes, tmp_long);
        
  lfBuf.lvalue = (lfBuf.lvalue >> ignBytesLength)<<ignBytesLength;
  
  //double tmpValue = lfBuf.value;
  
  vce->data = lfBuf.value+medianValue;
  vce->curValue = tmp_long;
  vce->reqBytesLength = reqBytesLength;
  vce->resiBitsLength = resiBitsLength;
}

size_t computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
  size_t dataLength;
  if(r1==0) 
  {
    dataLength = 0;
  }
  else if(r2==0) 
  {
    dataLength = r1;
  }
  else if(r3==0) 
  {
    dataLength = r1*r2;
  }
  else if(r4==0) 
  {
    dataLength = r1*r2*r3;
  }
  else if(r5==0) 
  {
    dataLength = r1*r2*r3*r4;
  }
  else 
  {
    dataLength = r1*r2*r3*r4*r5;
  }
  return dataLength;
}



float computeRangeSize_float(float* oriData, size_t size, float* valueRangeSize, float* medianValue)
{
  size_t i = 0;
  float min = oriData[0];
  float max = min;
  for(i=1;i<size;i++)
  {
    float data = oriData[i];
    if(min>data)
      min = data;
    else if(max<data)
      max = data;
  }

  *valueRangeSize = max - min;
  *medianValue = min + *valueRangeSize/2;
  return min;
}

double computeRangeSize_double(double* oriData, size_t size, double* valueRangeSize, double* medianValue)
{
  size_t i = 0;
  double min = oriData[0];
  double max = min;
  for(i=1;i<size;i++)
  {
    double data = oriData[i];
    if(min>data)
      min = data;
    else if(max<data)
      max = data;
  }
  
  *valueRangeSize = max - min;
  *medianValue = min + *valueRangeSize/2;
  return min;
}


long computeRangeSize_int(void* oriData, int dataType, size_t size, int64_t* valueRangeSize)
{
 size_t i = 0;
  long max = 0, min = 0;
  if(dataType == SZ_INT8)
  {
    char* data = (char*)oriData;
    char data_;
    min = data[0], max = min;
    computeMinMax(data);
  }
  else if(dataType == SZ_INT16)
  { 
    short* data = (short*)oriData;
    short data_; 
    min = data[0], max = min;
    computeMinMax(data);
  }
  else if(dataType == SZ_INT32)
  {
    int* data = (int*)oriData;
    int data_; 
    min = data[0], max = min;
    computeMinMax(data);
  }
  
  *valueRangeSize = max - min;
  return min; 
}


double getRealPrecision_float(float valueRangeSize, int errBoundMode, double absErrBound, double relBoundRatio, int *status)
{
  int state = SZ_SCES;
  double precision = 0;
  if(errBoundMode==REL)  
    precision = relBoundRatio*valueRangeSize;

  else
  {
   // printf("Error: error-bound-mode is incorrect!\n");
    state = SZ_BERR;
  }
  *status = state;
  return precision;
}


double getRealPrecision_int(long valueRangeSize, int errBoundMode, double absErrBound, double relBoundRatio, int *status)
{
  int state = SZ_SCES;
  double precision = 0;
 
  if(errBoundMode==REL)          
    precision = relBoundRatio*valueRangeSize;
    
  else
  {
    printf("Error: error-bound-mode is incorrect!\n");
    state = SZ_BERR;
  }
  *status = state;
  return precision;
}


double getRealPrecision_double(double valueRangeSize, int errBoundMode, double absErrBound, double relBoundRatio, int *status)
{
  int state = SZ_SCES;
  double precision = 0;
  
  if(errBoundMode==REL)  
    precision = relBoundRatio*valueRangeSize;

  else
  {
  //  printf("Error: error-bound-mode is incorrect!\n");
    state = SZ_BERR;
  }
  *status = state;
  return precision;
}



unsigned int roundUpToPowerOf2(unsigned int base)
{
  base -= 1;

  base = base | (base >> 1);
  base = base | (base >> 2);
  base = base | (base >> 4);
  base = base | (base >> 8);
  base = base | (base >> 16);

  return base + 1;
} 



void updateQuantizationInfo(int quant_intervals)
{
  exe_params->intvCapacity = quant_intervals;
  exe_params->intvRadius = quant_intervals/2;
} 



short getPrecisionReqLength_double(double precision)
{
  ldouble lbuf;
  lbuf.value = precision;
  long long lvalue = lbuf.lvalue;
  
  int expValue = (int)((lvalue & 0x7FF0000000000000) >> 52);
  expValue -= 1023;
//  unsigned char the1stManBit = (unsigned char)((lvalue & 0x0008000000000000) >> 51);
//  if(the1stManBit==1)
//    expValue--;
  return (short)expValue;
}

 


void computeReqLength_float(double realPrecision, short radExpo, int* reqLength, float* medianValue)
{
  short reqExpo = getPrecisionReqLength_double(realPrecision);
  *reqLength = 9+radExpo - reqExpo; //radExpo-reqExpo == reqMantiLength
  if(*reqLength<9)
    *reqLength = 9;
  if(*reqLength>32)
  { 
    *reqLength = 32;
    *medianValue = 0;
  }     
}


 void computeReqLength_double(double realPrecision, short radExpo, int* reqLength, double* medianValue)
{
  short reqExpo = getPrecisionReqLength_double(realPrecision);
  *reqLength = 12+radExpo - reqExpo; //radExpo-reqExpo == reqMantiLength
  if(*reqLength<12)
    *reqLength = 12;
  if(*reqLength>64)
  {
    *reqLength = 64;
    *medianValue = 0;
  }
}

unsigned int optimize_intervals_float_1D_opt(float *oriData, size_t dataLength, double realPrecision)
{ 
  size_t i = 0, radiusIndex;
  float pred_value = 0, pred_err;
  size_t *intervals = (size_t*)malloc(confparams_cpr->maxRangeRadius*sizeof(size_t));
  memset(intervals, 0, confparams_cpr->maxRangeRadius*sizeof(size_t));
  size_t totalSampleSize = 0;//dataLength/confparams_cpr->sampleDistance;

  float * data_pos = oriData + 2;
  while(data_pos - oriData < dataLength){
    totalSampleSize++;
    pred_value = data_pos[-1];
    pred_err = fabs(pred_value - *data_pos);
    radiusIndex = (unsigned long)((pred_err/realPrecision+1)/2);
    if(radiusIndex>=confparams_cpr->maxRangeRadius)
      radiusIndex = confparams_cpr->maxRangeRadius - 1;     
    intervals[radiusIndex]++;

    data_pos += confparams_cpr->sampleDistance;
  }
  //compute the appropriate number
  size_t targetCount = totalSampleSize*confparams_cpr->predThreshold;
  size_t sum = 0;
  for(i=0;i<confparams_cpr->maxRangeRadius;i++)
  {
    sum += intervals[i];
    if(sum>targetCount)
      break;
  }
  if(i>=confparams_cpr->maxRangeRadius)
    i = confparams_cpr->maxRangeRadius-1;
    
  unsigned int accIntervals = 2*(i+1);
  unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);
  
  if(powerOf2<32)
    powerOf2 = 32;
  
  free(intervals);
  return powerOf2;
}

unsigned int optimize_intervals_int8_1D(int8_t *oriData, size_t dataLength, double realPrecision)
{  
  size_t i = 0, radiusIndex;
  int64_t pred_value = 0, pred_err;
  size_t *intervals = (size_t*)malloc(confparams_cpr->maxRangeRadius*sizeof(size_t));
  memset(intervals, 0, confparams_cpr->maxRangeRadius*sizeof(size_t));
  size_t totalSampleSize = dataLength/confparams_cpr->sampleDistance;
  for(i=2;i<dataLength;i++)
  {
    if(i%confparams_cpr->sampleDistance==0)
    {
      //pred_value = 2*oriData[i-1] - oriData[i-2];
      pred_value = oriData[i-1];
      pred_err = llabs(pred_value - oriData[i]);
      radiusIndex = (uint64_t)((pred_err/realPrecision+1)/2);
      if(radiusIndex>=confparams_cpr->maxRangeRadius)
        radiusIndex = confparams_cpr->maxRangeRadius - 1;     
      intervals[radiusIndex]++;
    }
  }
  //compute the appropriate number
  size_t targetCount = totalSampleSize*confparams_cpr->predThreshold;
  size_t sum = 0;
  for(i=0;i<confparams_cpr->maxRangeRadius;i++)
  {
    sum += intervals[i];
    if(sum>targetCount)
      break;
  }
  if(i>=confparams_cpr->maxRangeRadius)
    i = confparams_cpr->maxRangeRadius-1;
    
  unsigned int accIntervals = 2*(i+1);
  unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);
  
  if(powerOf2<32)
    powerOf2 = 32;
  
  free(intervals);
  //printf("accIntervals=%d, powerOf2=%d\n", accIntervals, powerOf2);
  return powerOf2;
}

unsigned int optimize_intervals_int16_1D(int16_t *oriData, size_t dataLength, double realPrecision)
{  
  size_t i = 0, radiusIndex;
  int64_t pred_value = 0, pred_err;
  size_t *intervals = (size_t*)malloc(confparams_cpr->maxRangeRadius*sizeof(size_t));
  memset(intervals, 0, confparams_cpr->maxRangeRadius*sizeof(size_t));
  size_t totalSampleSize = dataLength/confparams_cpr->sampleDistance;
  for(i=2;i<dataLength;i++)
  {
    if(i%confparams_cpr->sampleDistance==0)
    {
      //pred_value = 2*oriData[i-1] - oriData[i-2];
      pred_value = oriData[i-1];
      pred_err = llabs(pred_value - oriData[i]);
      radiusIndex = (uint64_t)((pred_err/realPrecision+1)/2);
      if(radiusIndex>=confparams_cpr->maxRangeRadius)
        radiusIndex = confparams_cpr->maxRangeRadius - 1;     
      intervals[radiusIndex]++;
    }
  }
  //compute the appropriate number
  size_t targetCount = totalSampleSize*confparams_cpr->predThreshold;
  size_t sum = 0;
  for(i=0;i<confparams_cpr->maxRangeRadius;i++)
  {
    sum += intervals[i];
    if(sum>targetCount)
      break;
  }
  if(i>=confparams_cpr->maxRangeRadius)
    i = confparams_cpr->maxRangeRadius-1;
    
  unsigned int accIntervals = 2*(i+1);
  unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);
  
  if(powerOf2<32)
    powerOf2 = 32;
  
  free(intervals);
  //printf("accIntervals=%d, powerOf2=%d\n", accIntervals, powerOf2);
  return powerOf2;
}


 short getExponent_float(float value)
{
  //int ivalue = floatToBigEndianInt(value);
  lfloat lbuf;
  lbuf.value = value;
  int ivalue = lbuf.ivalue;
  
  int expValue = (ivalue & 0x7F800000) >> 23;
  expValue -= 127;
  return (short)expValue;
}

unsigned int optimize_intervals_int32_1D(int32_t *oriData, size_t dataLength, double realPrecision)
{  
  size_t i = 0, radiusIndex;
  int64_t pred_value = 0, pred_err;
  size_t *intervals = (size_t*)malloc(confparams_cpr->maxRangeRadius*sizeof(size_t));
  memset(intervals, 0, confparams_cpr->maxRangeRadius*sizeof(size_t));
  size_t totalSampleSize = dataLength/confparams_cpr->sampleDistance;
  for(i=2;i<dataLength;i++)
  {
    if(i%confparams_cpr->sampleDistance==0)
    {
      //pred_value = 2*oriData[i-1] - oriData[i-2];
      pred_value = oriData[i-1];
      pred_err = llabs(pred_value - oriData[i]);
      radiusIndex = (uint64_t)((pred_err/realPrecision+1)/2);
      if(radiusIndex>=confparams_cpr->maxRangeRadius)
        radiusIndex = confparams_cpr->maxRangeRadius - 1;     
      intervals[radiusIndex]++;
    }
  }
  //compute the appropriate number
  size_t targetCount = totalSampleSize*confparams_cpr->predThreshold;
  size_t sum = 0;
  for(i=0;i<confparams_cpr->maxRangeRadius;i++)
  {
    sum += intervals[i];
    if(sum>targetCount)
      break;
  }
  if(i>=confparams_cpr->maxRangeRadius)
    i = confparams_cpr->maxRangeRadius-1;
    
  unsigned int accIntervals = 2*(i+1);
  unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);
  
  if(powerOf2<32)
    powerOf2 = 32;
  
  free(intervals);
  //printf("accIntervals=%d, powerOf2=%d\n", accIntervals, powerOf2);
  return powerOf2;
}


unsigned int optimize_intervals_double_1D_opt(double *oriData, size_t dataLength, double realPrecision)
{  
  size_t i = 0, radiusIndex;
  double pred_value = 0, pred_err;
  size_t *intervals = (size_t*)malloc(confparams_cpr->maxRangeRadius*sizeof(size_t));
  memset(intervals, 0, confparams_cpr->maxRangeRadius*sizeof(size_t));
  size_t totalSampleSize = 0;

  double * data_pos = oriData + 2;
  while(data_pos - oriData < dataLength){
    totalSampleSize++;
    pred_value = data_pos[-1];
    pred_err = fabs(pred_value - *data_pos);
    radiusIndex = (unsigned long)((pred_err/realPrecision+1)/2);
    if(radiusIndex>=confparams_cpr->maxRangeRadius)
      radiusIndex = confparams_cpr->maxRangeRadius - 1;     
    intervals[radiusIndex]++;

    data_pos += confparams_cpr->sampleDistance;
  }
  //compute the appropriate number
  size_t targetCount = totalSampleSize*confparams_cpr->predThreshold;
  size_t sum = 0;
  for(i=0;i<confparams_cpr->maxRangeRadius;i++)
  {
    sum += intervals[i];
    if(sum>targetCount)
      break;
  }
  if(i>=confparams_cpr->maxRangeRadius)
    i = confparams_cpr->maxRangeRadius-1;
    
  unsigned int accIntervals = 2*(i+1);
  unsigned int powerOf2 = roundUpToPowerOf2(accIntervals);
  
  if(powerOf2<32)
    powerOf2 = 32;
  
  free(intervals);
  return powerOf2;
}

 
TightDataPointStorageF* SZ_compress_float_1D_MDQ(float *oriData, 
size_t dataLength, double realPrecision, float valueRangeSize, float medianValue_f)
{
  
  unsigned int quantization_intervals;
  if(exe_params->optQuantMode==1)
    quantization_intervals = optimize_intervals_float_1D_opt(oriData, dataLength, realPrecision);
  else
    quantization_intervals = exe_params->intvCapacity;
  updateQuantizationInfo(quantization_intervals); 

  size_t i;
  int reqLength;
  float medianValue = medianValue_f;
  short radExpo = getExponent_float(valueRangeSize/2);
  
    computeReqLength_float(realPrecision, radExpo, &reqLength, &medianValue); 

  int* type = (int*) malloc(dataLength*sizeof(int));
    
  float* spaceFillingValue = oriData; //
  
  DynamicIntArray *exactLeadNumArray;
  new_DIA(&exactLeadNumArray, DynArrayInitLen);
  
  DynamicByteArray *exactMidByteArray;
  new_DBA(&exactMidByteArray, DynArrayInitLen);
  
  DynamicIntArray *resiBitArray;
  new_DIA(&resiBitArray, DynArrayInitLen);
  
  unsigned char preDataBytes[4];
  intToBytes_bigEndian(preDataBytes, 0);
  
  int reqBytesLength = reqLength/8;
  int resiBitsLength = reqLength%8;
  float last3CmprsData[3] = {0};

  FloatValueCompressElement *vce = (FloatValueCompressElement*)malloc(sizeof(FloatValueCompressElement));
  LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));
        
  //add the first data  
  type[0] = 0;
  compressSingleFloatValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
  updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
  memcpy(preDataBytes,vce->curBytes,4);
  addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
  listAdd_float(last3CmprsData, vce->data);

    
  //add the second data
  type[1] = 0;
  compressSingleFloatValue(vce, spaceFillingValue[1], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
  updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
  memcpy(preDataBytes,vce->curBytes,4);
  addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
  listAdd_float(last3CmprsData, vce->data);

  int state;
  double checkRadius;
  float curData;
  float pred;
  float predAbsErr;
  checkRadius = (exe_params->intvCapacity-1)*realPrecision;
  double interval = 2*realPrecision;
  
  for(i=2;i<dataLength;i++)
  { 
    curData = spaceFillingValue[i];
    //pred = 2*last3CmprsData[0] - last3CmprsData[1];
    pred = last3CmprsData[0];
    predAbsErr = fabs(curData - pred);  
    if(predAbsErr<checkRadius)
    {
      state = (predAbsErr/realPrecision+1)/2;
      if(curData>=pred)
      {
        type[i] = exe_params->intvRadius+state;
        pred = pred + state*interval;
      }
      else //curData<pred
      {
        type[i] = exe_params->intvRadius-state;
        pred = pred - state*interval;
      }
        
      //double-check the prediction error in case of machine-epsilon impact 
      if(fabs(curData-pred)>realPrecision)
      { 
        type[i] = 0;        
        compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
        updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
        memcpy(preDataBytes,vce->curBytes,4);
        addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);    
        
        listAdd_float(last3CmprsData, vce->data); 

      }
      else
      {
        listAdd_float(last3CmprsData, pred);

      } 
      continue;
    }
    
    //unpredictable data processing   
    type[i] = 0;    
    compressSingleFloatValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
    updateLossyCompElement_Float(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
    memcpy(preDataBytes,vce->curBytes,4);
    addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);

    listAdd_float(last3CmprsData, vce->data);

    
  }//end of for
    
//  char* expSegmentsInBytes;
//  int expSegmentsInBytes_size = convertESCToBytes(esc, &expSegmentsInBytes);
  size_t exactDataNum = exactLeadNumArray->size;
  
  TightDataPointStorageF* tdps;
      
  new_TightDataPointStorageF(&tdps, dataLength, exactDataNum, 
      type, exactMidByteArray->array, exactMidByteArray->size,  
      exactLeadNumArray->array,  
      resiBitArray->array, resiBitArray->size, 
      resiBitsLength,
      realPrecision, medianValue, (char)reqLength, quantization_intervals, NULL, 0, 0);

/*
  printf("type\n");
  for(i=0;i<dataLength;i++)
    printf("%d ",type[i]);
  printf("\n");
*/


  //free memory
  free_DIA(exactLeadNumArray);
  free_DIA(resiBitArray);
  free(type); 
  free(vce);
  free(lce);  
  free(exactMidByteArray); //exactMidByteArray->array has been released in free_TightDataPointStorageF(tdps);
  
  return tdps;
}


 short getExponent_double(double value)
{
  //long lvalue = doubleToBigEndianLong(value);
  
  ldouble lbuf;
  lbuf.value = value;
  long long lvalue = lbuf.lvalue;
  
  int expValue = (int)((lvalue & 0x7FF0000000000000) >> 52);
  expValue -= 1023;
  return (short)expValue;
}



TightDataPointStorageD* SZ_compress_double_1D_MDQ(double *oriData, 
size_t dataLength, double realPrecision, double valueRangeSize, double medianValue_d)
{
  /*
#ifdef HAVE_TIMECMPR
  double* decData = NULL; 
  if(confparams_cpr->szMode == SZ_TEMPORAL_COMPRESSION)
    decData = (double*)(multisteps->hist_data);
#endif  
  */
  unsigned int quantization_intervals;
  if(exe_params->optQuantMode==1)
    quantization_intervals = optimize_intervals_double_1D_opt(oriData, dataLength, realPrecision);
  else
    quantization_intervals = exe_params->intvCapacity;
  updateQuantizationInfo(quantization_intervals); 
  int intvRadius = quantization_intervals/2;

  size_t i;
  int reqLength;
  double medianValue = medianValue_d;
  short radExpo = getExponent_double(valueRangeSize/2);

  computeReqLength_double(realPrecision, radExpo, &reqLength, &medianValue);  

  int* type = (int*) malloc(dataLength*sizeof(int));
    
  double* spaceFillingValue = oriData; //
  
  DynamicIntArray *exactLeadNumArray;
  new_DIA(&exactLeadNumArray, DynArrayInitLen);
  
  DynamicByteArray *exactMidByteArray;
  new_DBA(&exactMidByteArray, DynArrayInitLen);
  
  DynamicIntArray *resiBitArray;
  new_DIA(&resiBitArray, DynArrayInitLen);

  unsigned char preDataBytes[8];
  longToBytes_bigEndian(preDataBytes, 0);
  
  int reqBytesLength = reqLength/8;
  int resiBitsLength = reqLength%8;
  double last3CmprsData[3] = {0};

  DoubleValueCompressElement *vce = (DoubleValueCompressElement*)malloc(sizeof(DoubleValueCompressElement));
  LossyCompressionElement *lce = (LossyCompressionElement*)malloc(sizeof(LossyCompressionElement));     
        
  //add the first data  
  type[0] = 0;
  compressSingleDoubleValue(vce, spaceFillingValue[0], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
  updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
  memcpy(preDataBytes,vce->curBytes,8);
  addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
  listAdd_double(last3CmprsData, vce->data);
  /*
#ifdef HAVE_TIMECMPR  
  if(confparams_cpr->szMode == SZ_TEMPORAL_COMPRESSION)
    decData[0] = vce->data;
#endif    
    */
  //add the second data
  type[1] = 0;
  compressSingleDoubleValue(vce, spaceFillingValue[1], realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
  updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
  memcpy(preDataBytes,vce->curBytes,8);
  addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
  listAdd_double(last3CmprsData, vce->data);
  
  int state;
  double checkRadius;
  double curData;
  double pred ; 
  double predAbsErr;
  checkRadius = (quantization_intervals-1)*realPrecision;
  double interval = 2*realPrecision;

  double recip_realPrecision = 1/realPrecision;
  for(i=2;i<dataLength;i++)
  {       
    //printf("%.30G\n",last3CmprsData[0]);
    curData = spaceFillingValue[i];
    //pred = 2*last3CmprsData[0] - last3CmprsData[1];
    pred = last3CmprsData[0];
    predAbsErr = fabs(curData - pred);  
    if(predAbsErr<checkRadius)
    {
      state = (predAbsErr*recip_realPrecision+1)*0.5;
    //    state = (predAbsErr/realPrecision+1)/2;

      if(curData>=pred)
      {
        type[i] = intvRadius+state;
      //  type[i] = exe_params->intvRadius+state;

        pred = pred + state*interval;
      }
      else //curData<pred
      {
        type[i] = intvRadius-state;
       // type[i] = exe_params->intvRadius-state;

        pred = pred - state*interval;
      }
    //  listAdd_double(last3CmprsData, pred);
    /*
#ifdef HAVE_TIMECMPR                                                                    
      if(confparams_cpr->szMode == SZ_TEMPORAL_COMPRESSION)
        decData[i] = pred;      
#endif  
*/
//double-check the prediction error in case of machine-epsilon impact 
      if(fabs(curData-pred)>realPrecision)
      { 
        type[i] = 0;    
        compressSingleDoubleValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
        updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
        memcpy(preDataBytes,vce->curBytes,8);
        addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
              
        listAdd_double(last3CmprsData, vce->data); 

      }
      else
      {
        listAdd_double(last3CmprsData, pred);

      } 
      
      continue;
    }
    
    //unpredictable data processing
    type[i] = 0;    
    compressSingleDoubleValue(vce, curData, realPrecision, medianValue, reqLength, reqBytesLength, resiBitsLength);
    updateLossyCompElement_Double(vce->curBytes, preDataBytes, reqBytesLength, resiBitsLength, lce);
    memcpy(preDataBytes,vce->curBytes,8);
    addExactData(exactMidByteArray, exactLeadNumArray, resiBitArray, lce);
              
    listAdd_double(last3CmprsData, vce->data);
    /*
    pred = vce->data;
    
#ifdef HAVE_TIMECMPR
    if(confparams_cpr->szMode == SZ_TEMPORAL_COMPRESSION)
      decData[i] = vce->data;
#endif  
  */  
  }//end of for
    
  size_t exactDataNum = exactLeadNumArray->size;
  
  TightDataPointStorageD* tdps;
      
  new_TightDataPointStorageD(&tdps, dataLength, exactDataNum, 
      type, exactMidByteArray->array, exactMidByteArray->size,  
      exactLeadNumArray->array,  
      resiBitArray->array, resiBitArray->size, 
      resiBitsLength, 
      realPrecision, medianValue, (char)reqLength, quantization_intervals, NULL, 0, 0);
  
//  printf("exactDataNum=%d, expSegmentsInBytes_size=%d, exactMidByteArray->size=%d\n", 
//      exactDataNum, expSegmentsInBytes_size, exactMidByteArray->size);
  
  //free memory
  free_DIA(exactLeadNumArray);
  free_DIA(resiBitArray);
  free(type);
  free(vce);
  free(lce);  
  free(exactMidByteArray); //exactMidByteArray->array has been released in free_TightDataPointStorageF(tdps); 
  
  return tdps;  
}


void memcpyDBA_Data(DynamicByteArray *dba, unsigned char* data, size_t length)
{
  if(dba->size + length > dba->capacity)
  {
    dba->capacity = dba->size + length;
    dba->array = (unsigned char *)realloc(dba->array, dba->capacity*sizeof(unsigned char));
  }
  memcpy(&(dba->array[dba->size]), data, length);
  dba->size += length;
}

int computeByteSizePerIntValue(long valueRangeSize)
{
  if(valueRangeSize<=256)
    return 1;
  else if(valueRangeSize<=65536)
    return 2;
  else if(valueRangeSize<=4294967296) //2^32
    return 4;
  else
    return 8;
}


TightDataPointStorageI* SZ_compress_int8_1D_MDQ(int8_t *oriData, size_t dataLength, double realPrecision, int64_t valueRangeSize, int64_t minValue)
{
  unsigned char bytes[8] = {0,0,0,0,0,0,0,0};
  int byteSize = computeByteSizePerIntValue(valueRangeSize);
  unsigned int quantization_intervals;
  if(exe_params->optQuantMode==1)
    quantization_intervals = optimize_intervals_int8_1D(oriData, dataLength, realPrecision);
  else
    quantization_intervals = exe_params->intvCapacity;
  updateQuantizationInfo(quantization_intervals); 
  size_t i;

  int* type = (int*) malloc(dataLength*sizeof(int));
    
  int8_t* spaceFillingValue = oriData; //
  
  DynamicByteArray *exactDataByteArray;
  new_DBA(&exactDataByteArray, DynArrayInitLen);
    
  int64_t last3CmprsData[3] = {0,0,0};
        
  //add the first data  
  type[0] = 0;
  compressInt8Value(spaceFillingValue[0], minValue, byteSize, bytes);
  memcpyDBA_Data(exactDataByteArray, bytes, byteSize);
  listAdd_int(last3CmprsData, spaceFillingValue[0]);
    
  type[1] = 0;
  compressInt8Value(spaceFillingValue[1], minValue, byteSize, bytes);
  memcpyDBA_Data(exactDataByteArray, bytes, byteSize);
  listAdd_int(last3CmprsData, spaceFillingValue[1]);
  //printf("%.30G\n",last3CmprsData[0]);  
  
  int state;
  double checkRadius = (exe_params->intvCapacity-1)*realPrecision;
  int64_t curData;
  int64_t pred, predAbsErr;
  double interval = 2*realPrecision;
  
  for(i=2;i<dataLength;i++)
  {
    curData = spaceFillingValue[i];
    //pred = 2*last3CmprsData[0] - last3CmprsData[1];
    pred = last3CmprsData[0];
    predAbsErr = llabs(curData - pred); 
    if(predAbsErr<checkRadius)
    {
      state = (predAbsErr/realPrecision+1)/2;
      if(curData>=pred)
      {
        type[i] = exe_params->intvRadius+state;
        pred = pred + state*interval;
      }
      else //curData<pred
      {
        type[i] = exe_params->intvRadius-state;
        pred = pred - state*interval;
      }
      if(pred>SZ_INT8_MAX) pred = SZ_INT8_MAX;
      if(pred<SZ_INT8_MIN) pred = SZ_INT8_MIN;      
      listAdd_int(last3CmprsData, pred);          
      continue;
    }
    
    //unpredictable data processing   
    type[i] = 0;
    compressInt8Value(curData, minValue, byteSize, bytes);
    memcpyDBA_Data(exactDataByteArray, bytes, byteSize);
    listAdd_int(last3CmprsData, curData);
  }//end of for
    
  size_t exactDataNum = exactDataByteArray->size / byteSize;
  
  TightDataPointStorageI* tdps; 
      
  new_TightDataPointStorageI(&tdps, dataLength, exactDataNum, byteSize, 
      type, exactDataByteArray->array, exactDataByteArray->size,  
      realPrecision, minValue, quantization_intervals, SZ_INT8);

  //free memory
  free(type); 
  free(exactDataByteArray); 
  
  return tdps;
}



TightDataPointStorageI* SZ_compress_int16_1D_MDQ(int16_t *oriData, size_t dataLength, double realPrecision, int64_t valueRangeSize, int64_t minValue)
{
  unsigned char bytes[8] = {0,0,0,0,0,0,0,0};
  int byteSize = computeByteSizePerIntValue(valueRangeSize);
  unsigned int quantization_intervals;
  if(exe_params->optQuantMode==1)
    quantization_intervals = optimize_intervals_int16_1D(oriData, dataLength, realPrecision);
  else
    quantization_intervals = exe_params->intvCapacity;
  updateQuantizationInfo(quantization_intervals); 
  size_t i;

  int* type = (int*) malloc(dataLength*sizeof(int));
    
  int16_t* spaceFillingValue = oriData; //
  
  DynamicByteArray *exactDataByteArray;
  new_DBA(&exactDataByteArray, DynArrayInitLen);
    
  int64_t last3CmprsData[3] = {0,0,0};
        
  //add the first data  
  type[0] = 0;
  compressInt16Value(spaceFillingValue[0], minValue, byteSize, bytes);
  memcpyDBA_Data(exactDataByteArray, bytes, byteSize);
  listAdd_int(last3CmprsData, spaceFillingValue[0]);
    
  type[1] = 0;
  compressInt16Value(spaceFillingValue[1], minValue, byteSize, bytes);
  memcpyDBA_Data(exactDataByteArray, bytes, byteSize);
  listAdd_int(last3CmprsData, spaceFillingValue[1]);
  //printf("%.30G\n",last3CmprsData[0]);  
  
  int state;
  double checkRadius = (exe_params->intvCapacity-1)*realPrecision;
  int64_t curData;
  int64_t pred, predAbsErr;
  double interval = 2*realPrecision;
  
  for(i=2;i<dataLength;i++)
  {
    curData = spaceFillingValue[i];
    //pred = 2*last3CmprsData[0] - last3CmprsData[1];
    pred = last3CmprsData[0];
    predAbsErr = llabs(curData - pred); 
    if(predAbsErr<checkRadius)
    {
      state = (predAbsErr/realPrecision+1)/2;
      if(curData>=pred)
      {
        type[i] = exe_params->intvRadius+state;
        pred = pred + state*interval;
      }
      else //curData<pred
      {
        type[i] = exe_params->intvRadius-state;
        pred = pred - state*interval;
      }
      if(pred>SZ_INT16_MAX) pred = SZ_INT16_MAX;
      if(pred<SZ_INT16_MIN) pred = SZ_INT16_MIN;      
      listAdd_int(last3CmprsData, pred);          
      continue;
    }
    
    //unpredictable data processing   
    type[i] = 0;
    compressInt16Value(curData, minValue, byteSize, bytes);
    memcpyDBA_Data(exactDataByteArray, bytes, byteSize);
    listAdd_int(last3CmprsData, curData);
  }//end of for
    
  size_t exactDataNum = exactDataByteArray->size / byteSize;
  
  TightDataPointStorageI* tdps; 
      
  new_TightDataPointStorageI(&tdps, dataLength, exactDataNum, byteSize, 
      type, exactDataByteArray->array, exactDataByteArray->size,  
      realPrecision, minValue, quantization_intervals, SZ_INT16);

//sdi:Debug
/*  int sum =0;
  for(i=0;i<dataLength;i++)
    if(type[i]==0) sum++;
  printf("opt_quantizations=%d, exactDataNum=%d, sum=%d\n",quantization_intervals, exactDataNum, sum);*/
  
  //free memory
  free(type); 
  free(exactDataByteArray); //exactDataByteArray->array has been released in free_TightDataPointStorageF(tdps);
  
  return tdps;
}


TightDataPointStorageI* SZ_compress_int32_1D_MDQ(int32_t *oriData, size_t dataLength, double realPrecision, int64_t valueRangeSize, int64_t minValue)
{
  unsigned char bytes[8] = {0,0,0,0,0,0,0,0};
  int byteSize = computeByteSizePerIntValue(valueRangeSize);
  unsigned int quantization_intervals;
  if(exe_params->optQuantMode==1)
    quantization_intervals = optimize_intervals_int32_1D(oriData, dataLength, realPrecision);
  else
    quantization_intervals = exe_params->intvCapacity;
  updateQuantizationInfo(quantization_intervals); 
  size_t i;

  int* type = (int*) malloc(dataLength*sizeof(int));
    
  int32_t* spaceFillingValue = oriData; //
  
  DynamicByteArray *exactDataByteArray;
  new_DBA(&exactDataByteArray, DynArrayInitLen);
    
  int64_t last3CmprsData[3] = {0,0,0};
        
  //add the first data  
  type[0] = 0;
  compressInt32Value(spaceFillingValue[0], minValue, byteSize, bytes);
  memcpyDBA_Data(exactDataByteArray, bytes, byteSize);
  listAdd_int(last3CmprsData, spaceFillingValue[0]);
    
  type[1] = 0;
  compressInt32Value(spaceFillingValue[1], minValue, byteSize, bytes);
  memcpyDBA_Data(exactDataByteArray, bytes, byteSize);
  listAdd_int(last3CmprsData, spaceFillingValue[1]);
  //printf("%.30G\n",last3CmprsData[0]);  
  
  int state;
  double checkRadius = (exe_params->intvCapacity-1)*realPrecision;
  int64_t curData;
  int32_t pred, predAbsErr;
  double interval = 2*realPrecision;
  
  for(i=2;i<dataLength;i++)
  {

    curData = spaceFillingValue[i];
    pred = last3CmprsData[0];
    predAbsErr = llabs(curData - pred); 
    if(predAbsErr<checkRadius)
    {
      state = (predAbsErr/realPrecision+1)/2;
      if(curData>=pred)
      {
        type[i] = exe_params->intvRadius+state;
        pred = pred + state*interval;
      }
      else //curData<pred
      {
        type[i] = exe_params->intvRadius-state;
        pred = pred - state*interval;
      }

      listAdd_int(last3CmprsData, pred);          
      continue;
    }
    
    //unpredictable data processing   
    type[i] = 0;
    compressInt32Value(curData, minValue, byteSize, bytes);
    memcpyDBA_Data(exactDataByteArray, bytes, byteSize);
    listAdd_int(last3CmprsData, curData);
  }//end of for
    
  size_t exactDataNum = exactDataByteArray->size / byteSize;
  
  TightDataPointStorageI* tdps; 
      
  new_TightDataPointStorageI(&tdps, dataLength, exactDataNum, byteSize, 
      type, exactDataByteArray->array, exactDataByteArray->size,  
      realPrecision, minValue, quantization_intervals, SZ_INT32);

  //free memory
  free(type); 
  free(exactDataByteArray); //exactDataByteArray->array has been released in free_TightDataPointStorageF(tdps);
  
  return tdps;
}


void convertSZParamsToBytes(sz_params* params, unsigned char* result)
{
  unsigned char buf;
  buf = exe_params->optQuantMode;
  buf = (buf << 1) | dataEndianType;
  buf = (buf << 1) | sysEndianType;
  buf = (buf << 1) | params->szMode;
  
  int tmp = 0;
  switch(params->gzipMode)
  {
  case Z_BEST_SPEED:
    tmp = 0;
    break;
  case Z_DEFAULT_STRATEGY:
    tmp = 1;
    break;
  case Z_BEST_COMPRESSION:
    tmp = 2;
    break;
    }
  buf = (buf << 2) | tmp;
  buf = (buf << 2) |  params->pwr_type;
  result[0] = buf;
  
    //sampleDistance; //2 bytes
    int16ToBytes_bigEndian(&result[1], params->sampleDistance);
    
    //conf_params->predThreshold;  // 2 bytes
    short tmp2 = params->predThreshold * 10000;
    int16ToBytes_bigEndian(&result[3], tmp2);
     
    //errorBoundMode; //4bits(0.5 byte)
    result[5] = params->errorBoundMode;
    result[5] = (result[5] << 4) | (params->dataType & 0x17);
     
    switch(params->errorBoundMode)
    {

  case REL:
    memset(&result[6], 0, 4);
    floatToBytes(&result[10], (float)(params->relBoundRatio)); //big_endian
    break;
 
  }
   
    //segment_size  // 2 bytes
    int16ToBytes_bigEndian(&result[14], (short)(params->segment_size));
    
    if(exe_params->optQuantMode==1)
    int32ToBytes_bigEndian(&result[16], params->max_quant_intervals);
  else
    int32ToBytes_bigEndian(&result[16], params->quantization_intervals);
}



void convertTDPStoBytes_int(TightDataPointStorageI* tdps, unsigned char* bytes, unsigned char sameByte)
{
  size_t i, k = 0;
  
  unsigned char byteBuffer[8] = {0,0,0,0,0,0,0,0};
  
  for(i = 0;i<3;i++)//3 bytes
    bytes[k++] = 0;  //versionNumber[i];
  bytes[k++] = sameByte;  //1 byte
  
  convertSZParamsToBytes(confparams_cpr, &(bytes[k]));
  k = k + MetaDataByteLength; 
    
  bytes[k++] = tdps->exactByteSize; //1 byte

  sizeToBytes(byteBuffer, tdps->dataSeriesLength);
  for(i = 0;i<exe_params->SZ_SIZE_TYPE;i++)//ST: 4 or 8 bytes
    bytes[k++] = byteBuffer[i]; 
  
  intToBytes_bigEndian(byteBuffer, confparams_cpr->max_quant_intervals);
  for(i = 0;i<4;i++)//4
    bytes[k++] = byteBuffer[i];
  
  intToBytes_bigEndian(byteBuffer, tdps->intervals);
  for(i = 0;i<4;i++)//4
    bytes[k++] = byteBuffer[i];     
  
  longToBytes_bigEndian(byteBuffer, tdps->minValue);
  for (i = 0; i < 8; i++)// 8
    bytes[k++] = byteBuffer[i];

  doubleToBytes(byteBuffer, tdps->realPrecision);
  for (i = 0; i < 8; i++)// 8
    bytes[k++] = byteBuffer[i];     

  sizeToBytes(byteBuffer, tdps->typeArray_size);
  for(i = 0;i<exe_params->SZ_SIZE_TYPE;i++)//ST
    bytes[k++] = byteBuffer[i];

  sizeToBytes(byteBuffer, tdps->exactDataNum);
  for(i = 0;i<exe_params->SZ_SIZE_TYPE;i++)//ST
    bytes[k++] = byteBuffer[i];

  sizeToBytes(byteBuffer, tdps->exactDataBytes_size);
  for(i = 0;i<exe_params->SZ_SIZE_TYPE;i++)//ST
    bytes[k++] = byteBuffer[i];

  memcpy(&(bytes[k]), tdps->typeArray, tdps->typeArray_size);
  k += tdps->typeArray_size;

  memcpy(&(bytes[k]), tdps->exactDataBytes, tdps->exactDataBytes_size);
  k += tdps->exactDataBytes_size;
}

void convertTDPStoBytes_float(TightDataPointStorageF* tdps, unsigned char* bytes, unsigned char* dsLengthBytes, unsigned char sameByte)
{
  size_t i, k = 0;
  unsigned char intervalsBytes[4];
  unsigned char typeArrayLengthBytes[8];
  unsigned char exactLengthBytes[8];
  unsigned char exactMidBytesLength[8];
  unsigned char realPrecisionBytes[8];
  
  unsigned char medianValueBytes[4];
  
  unsigned char segment_sizeBytes[8];
  unsigned char pwrErrBoundBytes_sizeBytes[4];
  unsigned char max_quant_intervals_Bytes[4];
  
  
  for(i = 0;i<3;i++)//3 bytes
    bytes[k++] = 0;//versionNumber[i];
  bytes[k++] = sameByte;  //1 byte  
  
  convertSZParamsToBytes(confparams_cpr, &(bytes[k]));
//  k = k +MetaDataByteLength;
  
  for(i = 0;i<exe_params->SZ_SIZE_TYPE;i++)//ST: 4 or 8 bytes
    bytes[k++] = dsLengthBytes[i];  
  intToBytes_bigEndian(max_quant_intervals_Bytes, confparams_cpr->max_quant_intervals);
  for(i = 0;i<4;i++)//4
    bytes[k++] = max_quant_intervals_Bytes[i];    
  
/*  if(confparams_cpr->errorBoundMode>=PW_REL)
  {
    bytes[k++] = tdps->radExpo; //1 byte      
    
    sizeToBytes(segment_sizeBytes, confparams_cpr->segment_size);
    for(i = 0;i<exe_params->SZ_SIZE_TYPE;i++)//ST
      bytes[k++] = segment_sizeBytes[i];        
      
    intToBytes_bigEndian(pwrErrBoundBytes_sizeBytes, tdps->pwrErrBoundBytes_size);
    for(i = 0;i<4;i++)//4
      bytes[k++] = pwrErrBoundBytes_sizeBytes[i];         
  }
*/  
  intToBytes_bigEndian(intervalsBytes, tdps->intervals);
  for(i = 0;i<4;i++)//4
    bytes[k++] = intervalsBytes[i];     
  
  floatToBytes(medianValueBytes, tdps->medianValue);
  for (i = 0; i < 4; i++)// 4
    bytes[k++] = medianValueBytes[i];   

  bytes[k++] = tdps->reqLength; //1 byte

/*  if(errorBoundMode>=PW_REL)
    doubleToBytes(realPrecisionBytes, pw_relBoundRatio);
  else*/


  doubleToBytes(realPrecisionBytes, tdps->realPrecision);

  for (i = 0; i < 8; i++)// 8
    bytes[k++] = realPrecisionBytes[i];




    sizeToBytes(typeArrayLengthBytes, tdps->typeArray_size);
  for(i = 0;i<exe_params->SZ_SIZE_TYPE;i++)//ST
    bytes[k++] = typeArrayLengthBytes[i];




    sizeToBytes(exactLengthBytes, tdps->exactDataNum);
  for(i = 0;i<exe_params->SZ_SIZE_TYPE;i++)//ST
    bytes[k++] = exactLengthBytes[i];

  sizeToBytes(exactMidBytesLength, tdps->exactMidBytes_size);
  for(i = 0;i<exe_params->SZ_SIZE_TYPE;i++)//ST
    bytes[k++] = exactMidBytesLength[i];


  memcpy(&(bytes[k]), tdps->typeArray, tdps->typeArray_size);
  k += tdps->typeArray_size;

  memcpy(&(bytes[k]), tdps->leadNumArray, tdps->leadNumArray_size);
  k += tdps->leadNumArray_size;
  memcpy(&(bytes[k]), tdps->exactMidBytes, tdps->exactMidBytes_size);
  k += tdps->exactMidBytes_size;

  if(tdps->residualMidBits!=NULL)
  {
    memcpy(&(bytes[k]), tdps->residualMidBits, tdps->residualMidBits_size);
    k += tdps->residualMidBits_size;
  } 
}

void convertTDPStoBytes_double(TightDataPointStorageD* tdps, unsigned char* bytes, unsigned char* dsLengthBytes, unsigned char sameByte)
{
  size_t i, k = 0;
  unsigned char intervalsBytes[4];
  unsigned char typeArrayLengthBytes[8];
  unsigned char exactLengthBytes[8];
  unsigned char exactMidBytesLength[8];
  unsigned char realPrecisionBytes[8];
  
  unsigned char medianValueBytes[8];
  
  unsigned char segment_sizeBytes[8];
  unsigned char pwrErrBoundBytes_sizeBytes[4];
  unsigned char max_quant_intervals_Bytes[4];
  
  for(i = 0;i<3;i++)//3 bytes
    bytes[k++] = 0;  // versionNumber[i];
  bytes[k++] = sameByte;  //1 byte  
  
  convertSZParamsToBytes(confparams_cpr, &(bytes[k]));
  k = k + MetaDataByteLength_double;
  
  for(i = 0;i<exe_params->SZ_SIZE_TYPE;i++)//ST: 4 or 8 bytes
    bytes[k++] = dsLengthBytes[i];  
  intToBytes_bigEndian(max_quant_intervals_Bytes, confparams_cpr->max_quant_intervals);
  for(i = 0;i<4;i++)//4
    bytes[k++] = max_quant_intervals_Bytes[i];    
  
  intToBytes_bigEndian(intervalsBytes, tdps->intervals);
  for(i = 0;i<4;i++)//4
    bytes[k++] = intervalsBytes[i];   
  
  doubleToBytes(medianValueBytes, tdps->medianValue);
  for (i = 0; i < 8; i++)// 8
    bytes[k++] = medianValueBytes[i];   

  bytes[k++] = tdps->reqLength; //1 byte

  doubleToBytes(realPrecisionBytes, tdps->realPrecision);
  for (i = 0; i < 8; i++)// 8
    bytes[k++] = realPrecisionBytes[i];
      
  sizeToBytes(typeArrayLengthBytes, tdps->typeArray_size);
  for(i = 0;i<exe_params->SZ_SIZE_TYPE;i++)//ST
    bytes[k++] = typeArrayLengthBytes[i];       
        
  sizeToBytes(exactLengthBytes, tdps->exactDataNum);
  for(i = 0;i<exe_params->SZ_SIZE_TYPE;i++)//ST
    bytes[k++] = exactLengthBytes[i];

  sizeToBytes(exactMidBytesLength, tdps->exactMidBytes_size);
  for(i = 0;i<exe_params->SZ_SIZE_TYPE;i++)//ST
    bytes[k++] = exactMidBytesLength[i];
  
  memcpy(&(bytes[k]), tdps->typeArray, tdps->typeArray_size);
  k += tdps->typeArray_size;

  memcpy(&(bytes[k]), tdps->leadNumArray, tdps->leadNumArray_size);
  k += tdps->leadNumArray_size;
  memcpy(&(bytes[k]), tdps->exactMidBytes, tdps->exactMidBytes_size);
  k += tdps->exactMidBytes_size;

  if(tdps->residualMidBits!=NULL)
  {
    memcpy(&(bytes[k]), tdps->residualMidBits, tdps->residualMidBits_size);
    k += tdps->residualMidBits_size;
  }   
}

int convertDataTypeSizeCode(int dataTypeSizeCode)
{
  int result = 0;
  switch(dataTypeSizeCode)
  {
  case 0:
    result = 1;
    break;
  case 1:
    result = 2;
    break;
  case 2:
    result = 4;
    break;
  case 3:
    result = 8;
    break;
  }
  return result;  
}

int convertDataTypeSize(int dataTypeSize)
{
  int result = 0;
  switch(dataTypeSize)
  {
  case 1:
    result = 0; //0000
    break;
  case 2:
    result = 4; //0100
    break;
  case 4:
    result = 8; //1000
    break;
  case 8:
    result = 12; //1100
    break;
  }
  return result;
}



void convertTDPStoFlatBytes_float(TightDataPointStorageF *tdps, unsigned char** bytes, size_t *size)
{
  size_t i, k = 0; 
  unsigned char dsLengthBytes[8];
  
  if(exe_params->SZ_SIZE_TYPE==4)
    intToBytes_bigEndian(dsLengthBytes, tdps->dataSeriesLength);//4
  else
    longToBytes_bigEndian(dsLengthBytes, tdps->dataSeriesLength);//8
    
  unsigned char sameByte = tdps->allSameData==1?(unsigned char)1:(unsigned char)0;
  sameByte = sameByte | (confparams_cpr->szMode << 1);
  if(tdps->isLossless)
    sameByte = (unsigned char) (sameByte | 0x10);

  if(exe_params->SZ_SIZE_TYPE==8)
    sameByte = (unsigned char) (sameByte | 0x40); // 01000000, the 6th bit


if (tdps->rtypeArray == NULL)
  {
    size_t residualMidBitsLength = tdps->residualMidBits == NULL ? 0 : tdps->residualMidBits_size;
    size_t segmentL = 0, radExpoL = 0, pwrBoundArrayL = 0;
    int minLogValueSize = 0;

    size_t totalByteLength = 3 + 1 + /*MetaDataByteLength +*/ exe_params->SZ_SIZE_TYPE + 4 + radExpoL + segmentL + pwrBoundArrayL + 4 + 4 + 1 + 8 
        + exe_params->SZ_SIZE_TYPE + exe_params->SZ_SIZE_TYPE + exe_params->SZ_SIZE_TYPE + minLogValueSize
        + tdps->typeArray_size + tdps->leadNumArray_size 
        + tdps->exactMidBytes_size + residualMidBitsLength + tdps->pwrErrBoundBytes_size;

    //printf("totalByteLenght %d\n",totalByteLength);

    *bytes = (unsigned char *)malloc(sizeof(unsigned char)*totalByteLength);

    convertTDPStoBytes_float(tdps, *bytes, dsLengthBytes, sameByte);
    
    *size = totalByteLength;
  }

}

//convert TightDataPointStorageI to bytes...

void convertTDPStoFlatBytes_int(TightDataPointStorageI *tdps, unsigned char** bytes, size_t *size)
{
  size_t i, k = 0; 
  unsigned char dsLengthBytes[8];
  
  if(exe_params->SZ_SIZE_TYPE==4)
    intToBytes_bigEndian(dsLengthBytes, tdps->dataSeriesLength);//4
  else
    longToBytes_bigEndian(dsLengthBytes, tdps->dataSeriesLength);//8

  unsigned char sameByte = tdps->allSameData==1?(unsigned char)1:(unsigned char)0;
  sameByte = sameByte | (confparams_cpr->szMode << 1);
  if(tdps->isLossless)
    sameByte = (unsigned char) (sameByte | 0x10);
  
  if(exe_params->SZ_SIZE_TYPE==8)
    sameByte = (unsigned char) (sameByte | 0x40); // 01000000, the 6th bit
  
    if(confparams_cpr->errorBoundMode>=PW_REL)
    {     
      printf("Error: errorBoundMode >= PW_REL!! can't be...\n");
      exit(0);
    }

    size_t totalByteLength = 3 + 1 + MetaDataByteLength + 1 + exe_params->SZ_SIZE_TYPE + 4 + 4 + 8 + 8
        + exe_params->SZ_SIZE_TYPE + exe_params->SZ_SIZE_TYPE + exe_params->SZ_SIZE_TYPE
        + tdps->typeArray_size + tdps->exactDataBytes_size;

    *bytes = (unsigned char *)malloc(sizeof(unsigned char)*totalByteLength);
    memset(*bytes, 0, sizeof(unsigned char)*totalByteLength);

    convertTDPStoBytes_int(tdps, *bytes, sameByte);
    
    *size = totalByteLength;
  
}


//Convert TightDataPointStorageD to bytes...
void convertTDPStoFlatBytes_double(TightDataPointStorageD *tdps, unsigned char** bytes, size_t *size) 
{
  size_t i, k = 0; 
  unsigned char dsLengthBytes[8];
  
  if(exe_params->SZ_SIZE_TYPE==4)
    intToBytes_bigEndian(dsLengthBytes, tdps->dataSeriesLength);//4
  else
    longToBytes_bigEndian(dsLengthBytes, tdps->dataSeriesLength);//8
  
  unsigned char sameByte = tdps->allSameData==1?(unsigned char)1:(unsigned char)0;
  sameByte = sameByte | (confparams_cpr->szMode << 1);
  if(tdps->isLossless)
    sameByte = (unsigned char) (sameByte | 0x10); 
 
  if(exe_params->SZ_SIZE_TYPE==8)
    sameByte = (unsigned char) (sameByte | 0x40); // 01000000, the 6th bit
    
  if (tdps->rtypeArray == NULL) 
  {
    size_t residualMidBitsLength = tdps->residualMidBits == NULL ? 0 : tdps->residualMidBits_size;
    size_t segmentL = 0, radExpoL = 0, pwrBoundArrayL = 0;
    int minLogValueSize = 0;

    size_t totalByteLength = 3 + 1 + MetaDataByteLength_double + exe_params->SZ_SIZE_TYPE + 4 + radExpoL + segmentL + pwrBoundArrayL + 4 + 8 + 1 + 8 
        + exe_params->SZ_SIZE_TYPE + exe_params->SZ_SIZE_TYPE + exe_params->SZ_SIZE_TYPE 
        + minLogValueSize /*max absolute log value*/
        + tdps->typeArray_size + tdps->leadNumArray_size
        + tdps->exactMidBytes_size + residualMidBitsLength + tdps->pwrErrBoundBytes_size;
  
    *bytes = (unsigned char *)malloc(sizeof(unsigned char)*totalByteLength);

    convertTDPStoBytes_double(tdps, *bytes, dsLengthBytes, sameByte);
    
    *size = totalByteLength;
  }
  
}


 
char SZ_compress_args_float_NoCkRngeNoGzip_1D(unsigned char** newByteData, float *oriData, 
size_t dataLength, double realPrecision, size_t *outSize, float valueRangeSize, float medianValue_f)
{   
  char compressionType = 0; 
  TightDataPointStorageF* tdps = NULL;  

  tdps = SZ_compress_float_1D_MDQ(oriData, dataLength, realPrecision, valueRangeSize, medianValue_f); 

  convertTDPStoFlatBytes_float(tdps, newByteData, outSize);

  free_TightDataPointStorageF(tdps);
  return compressionType;
}


char SZ_compress_args_double_NoCkRngeNoGzip_1D(int cmprType, unsigned char** newByteData, double *oriData, 
size_t dataLength, double realPrecision, size_t *outSize, double valueRangeSize, double medianValue_d)
{
  char compressionType = 0; 
  TightDataPointStorageD* tdps = NULL;  
  tdps = SZ_compress_double_1D_MDQ(oriData, dataLength, realPrecision, valueRangeSize, medianValue_d);      
  convertTDPStoFlatBytes_double(tdps, newByteData, outSize);
  free_TightDataPointStorageD(tdps);  
  return compressionType;
}


void SZ_compress_args_int8_NoCkRngeNoGzip_1D(unsigned char** newByteData, int8_t *oriData, 
size_t dataLength, double realPrecision, size_t *outSize, int64_t valueRangeSize, int8_t minValue)
{
  TightDataPointStorageI* tdps = SZ_compress_int8_1D_MDQ(oriData, dataLength, realPrecision, valueRangeSize, minValue);
  //TODO: return bytes....
  convertTDPStoFlatBytes_int(tdps, newByteData, outSize);
 
  free_TightDataPointStorageI(tdps);
}

void SZ_compress_args_int16_NoCkRngeNoGzip_1D(unsigned char** newByteData, int16_t *oriData, 
size_t dataLength, double realPrecision, size_t *outSize, int64_t valueRangeSize, int16_t minValue)
{
  TightDataPointStorageI* tdps = SZ_compress_int16_1D_MDQ(oriData, dataLength, realPrecision, valueRangeSize, minValue);
  //TODO: return bytes....
  convertTDPStoFlatBytes_int(tdps, newByteData, outSize);
 
  free_TightDataPointStorageI(tdps);
}



void SZ_compress_args_int32_NoCkRngeNoGzip_1D(unsigned char** newByteData, int32_t *oriData, 
size_t dataLength, double realPrecision, size_t *outSize, int64_t valueRangeSize, int32_t minValue)
{
  TightDataPointStorageI* tdps = SZ_compress_int32_1D_MDQ(oriData, dataLength, realPrecision, valueRangeSize, minValue);
  //TODO: return bytes....
  convertTDPStoFlatBytes_int(tdps, newByteData, outSize);
  free_TightDataPointStorageI(tdps);
}

 
int SZ_compress_args_float(unsigned char** newByteData, float *oriData, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio, double pwRelBoundRatio)
{
  confparams_cpr->errorBoundMode = errBoundMode;

  int status = SZ_SCES;
  size_t dataLength = computeDataLength(r5,r4,r3,r2,r1);
  
  if(dataLength <= MIN_NUM_OF_ELEMENTS)
  {
//    printf("A VOIR\n");
    return status;
  }
  
  float valueRangeSize = 0, medianValue = 0;
  
  float min = computeRangeSize_float(oriData, dataLength, &valueRangeSize, &medianValue);
  float max = min+valueRangeSize;
  double realPrecision = 0; 
  

  realPrecision = getRealPrecision_float(valueRangeSize, errBoundMode, absErr_Bound, relBoundRatio, &status);
    
  {
    size_t tmpOutSize = 0;
    unsigned char* tmpByteData;
    
    if (r2==0)
    {

      {
        SZ_compress_args_float_NoCkRngeNoGzip_1D(&tmpByteData, oriData, r1, realPrecision, &tmpOutSize, valueRangeSize, medianValue);

       //     printf("tmpByteData size %d\n",tmpOutSize);
      }
    }

    else
    {
    //  printf("Error: doesn't support 5 dimensions for now.\n");
      status = SZ_DERR; //dimension error
    }
    //Call Gzip to do the further compression.
    if(confparams_cpr->szMode==SZ_BEST_SPEED)
    {
      *outSize = tmpOutSize;
      *newByteData = tmpByteData;
    }
    else if(confparams_cpr->szMode==SZ_BEST_COMPRESSION || confparams_cpr->szMode==SZ_DEFAULT_COMPRESSION || confparams_cpr->szMode==SZ_TEMPORAL_COMPRESSION)
    {
      *outSize = sz_lossless_compress(confparams_cpr->losslessCompressor, confparams_cpr->gzipMode, tmpByteData, tmpOutSize, newByteData);
            //printf("COMPRESSED SIZE %u\n",*outSize);
      free(tmpByteData);
    }
    else
    {
    //  printf("Error: Wrong setting of confparams_cpr->szMode in the float compression.\n");
      status = SZ_MERR; //mode error      
    }
  }
  
  return status;
}


double computeABSErrBoundFromPSNR(double psnr, double threshold, double value_range)
{
  double v1 = psnr + 10 * log10(1-2.0/3.0*threshold);
  double v2 = v1/(-20);
  double v3 = pow(10, v2);
  return value_range * v3;
} 

int SZ_compress_args_double(int cmprType, int withRegression, unsigned char** newByteData, double *oriData, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio, double pwRelBoundRatio)
{
 //confparams_cpr->dataType = SZ_DOUBLE;
  confparams_cpr->errorBoundMode = errBoundMode;
  
  int status = SZ_SCES;
  size_t dataLength = computeDataLength(r5,r4,r3,r2,r1);
  
  if(dataLength <= MIN_NUM_OF_ELEMENTS)
  {
   // *newByteData = SZ_skip_compress_double(oriData, dataLength, outSize);
    return status;
  }
  
  double valueRangeSize = 0, medianValue = 0;
  
  double  min = computeRangeSize_double(oriData, dataLength, &valueRangeSize, &medianValue);  
  double max = min+valueRangeSize;
  double realPrecision = 0;   
  realPrecision = getRealPrecision_double(valueRangeSize, errBoundMode, absErr_Bound, relBoundRatio, &status);
  {
    size_t tmpOutSize = 0;
    unsigned char* tmpByteData;
    if (r2==0)
    {      
            SZ_compress_args_double_NoCkRngeNoGzip_1D(cmprType, &tmpByteData, oriData, r1, realPrecision, &tmpOutSize, valueRangeSize, medianValue);
          //  printf("tmpByteData size %d\n",tmpOutSize);

    }
  
    else
    {
     // printf("Error: doesn't support 5 dimensions for now.\n");
      status = SZ_DERR;
    }
        
    //Call Gzip to do the further compression.
    if(confparams_cpr->szMode==SZ_BEST_SPEED)
    {
      *outSize = tmpOutSize;
      *newByteData = tmpByteData;     
    }
    else if(confparams_cpr->szMode==SZ_BEST_COMPRESSION || confparams_cpr->szMode==SZ_DEFAULT_COMPRESSION || confparams_cpr->szMode==SZ_TEMPORAL_COMPRESSION)
    {
      *outSize = sz_lossless_compress(confparams_cpr->losslessCompressor, confparams_cpr->gzipMode, tmpByteData, tmpOutSize, newByteData);
      free(tmpByteData);
    }
    else
    {
    //  printf("Error: Wrong setting of confparams_cpr->szMode in the double compression.\n");
      status = SZ_MERR; 
    }
  }

  return status;
}




int SZ_compress_args_int8(unsigned char** newByteData, int8_t *oriData, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio)
{
  confparams_cpr->errorBoundMode = errBoundMode;
  
  if(errBoundMode>=PW_REL)
  {
   // printf("Error: Current SZ version doesn't support integer data compression with point-wise relative error bound being based on pwrType=AVG\n");
    exit(0);
    return SZ_NSCS;
  }
  int status = SZ_SCES;
  size_t dataLength = computeDataLength(r5,r4,r3,r2,r1);
  int64_t valueRangeSize = 0;

  int8_t minValue = (int8_t)computeRangeSize_int(oriData, SZ_INT8, dataLength, &valueRangeSize);
  double realPrecision = 0; 
  
  if(confparams_cpr->errorBoundMode==PSNR)
  {
    confparams_cpr->errorBoundMode = ABS;
    realPrecision = confparams_cpr->absErrBound = computeABSErrBoundFromPSNR(confparams_cpr->psnr, (double)confparams_cpr->predThreshold, (double)valueRangeSize);
    //printf("realPrecision=%lf\n", realPrecision);
  }
  else
    realPrecision = getRealPrecision_int(valueRangeSize, errBoundMode, absErr_Bound, relBoundRatio, &status);
  
  {
    size_t tmpOutSize = 0;
    unsigned char* tmpByteData;
    if (r2==0)
    {
      SZ_compress_args_int8_NoCkRngeNoGzip_1D(&tmpByteData, oriData, r1, realPrecision, &tmpOutSize, valueRangeSize, minValue);
    }
    else
    {
    //  printf("Error: doesn't support 5 dimensions for now.\n");
      status = SZ_DERR; //dimension error
    }
    //Call Gzip to do the further compression.
    if(confparams_cpr->szMode==SZ_BEST_SPEED)
    {
      *outSize = tmpOutSize;
      *newByteData = tmpByteData;
    }
    else if(confparams_cpr->szMode==SZ_BEST_COMPRESSION || confparams_cpr->szMode==SZ_DEFAULT_COMPRESSION)
    {
      *outSize = sz_lossless_compress(confparams_cpr->losslessCompressor, confparams_cpr->gzipMode, tmpByteData, tmpOutSize, newByteData);
      free(tmpByteData);
    }
    else
    {
    //  printf("Error: Wrong setting of confparams_cpr->szMode in the int8_t compression.\n");
      status = SZ_MERR; //mode error      
    }
  }

  return status;
}

int SZ_compress_args_int16(unsigned char** newByteData, int16_t *oriData, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio)
{
  confparams_cpr->errorBoundMode = errBoundMode;
  
  if(errBoundMode>=PW_REL)
  {
  //  printf("Error: Current SZ version doesn't support integer data compression with point-wise relative error bound being based on pwrType=AVG\n");
    exit(0);
    return SZ_NSCS;
  }
  int status = SZ_SCES;
  size_t dataLength = computeDataLength(r5,r4,r3,r2,r1);
  int64_t valueRangeSize = 0;

  int16_t minValue = (int16_t)computeRangeSize_int(oriData, SZ_INT16, dataLength, &valueRangeSize);
  double realPrecision = 0; 
  
  if(confparams_cpr->errorBoundMode==PSNR)
  {
    confparams_cpr->errorBoundMode = ABS;
    realPrecision = confparams_cpr->absErrBound = computeABSErrBoundFromPSNR(confparams_cpr->psnr, (double)confparams_cpr->predThreshold, (double)valueRangeSize);
    //printf("realPrecision=%lf\n", realPrecision);
  }
  else
    realPrecision = getRealPrecision_int(valueRangeSize, errBoundMode, absErr_Bound, relBoundRatio, &status);

  {
    size_t tmpOutSize = 0;
    unsigned char* tmpByteData;
    if (r2==0)
    {
      SZ_compress_args_int16_NoCkRngeNoGzip_1D(&tmpByteData, oriData, r1, realPrecision, &tmpOutSize, valueRangeSize, minValue);
    }
    
    else
    {
    //  printf("Error: doesn't support 5 dimensions for now.\n");
      status = SZ_DERR; //dimension error
    }
    //Call Gzip to do the further compression.
    if(confparams_cpr->szMode==SZ_BEST_SPEED)
    {
      *outSize = tmpOutSize;
      *newByteData = tmpByteData;
    }
    else if(confparams_cpr->szMode==SZ_BEST_COMPRESSION || confparams_cpr->szMode==SZ_DEFAULT_COMPRESSION)
    {
      *outSize = sz_lossless_compress(confparams_cpr->losslessCompressor, confparams_cpr->gzipMode, tmpByteData, tmpOutSize, newByteData);
      free(tmpByteData);
    }
    else
    {
    //  printf("Error: Wrong setting of confparams_cpr->szMode in the int16_t compression.\n");
      status = SZ_MERR; //mode error      
    }
  }
  
  return status;
}


int SZ_compress_args_int32(unsigned char** newByteData, int32_t *oriData, 
size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t *outSize, 
int errBoundMode, double absErr_Bound, double relBoundRatio)
{
  confparams_cpr->errorBoundMode = errBoundMode;
  
  if(errBoundMode>=PW_REL)
  {
  //  printf("Error: Current SZ version doesn't support integer data compression with point-wise relative error bound being based on pwrType=AVG\n");
    exit(0);
    return SZ_NSCS;
  }
  int status = SZ_SCES;
  size_t dataLength = computeDataLength(r5,r4,r3,r2,r1);
  int64_t valueRangeSize = 0;

  int32_t minValue = (int32_t)computeRangeSize_int(oriData, SZ_INT32, dataLength, &valueRangeSize);
  double realPrecision = 0; 
  
  if(confparams_cpr->errorBoundMode==PSNR)
  {
    confparams_cpr->errorBoundMode = ABS;
    realPrecision = confparams_cpr->absErrBound = computeABSErrBoundFromPSNR(confparams_cpr->psnr, (double)confparams_cpr->predThreshold, (double)valueRangeSize);
    //printf("realPrecision=%lf\n", realPrecision);
  }
  else
    realPrecision = getRealPrecision_int(valueRangeSize, errBoundMode, absErr_Bound, relBoundRatio, &status);

  { 
    size_t tmpOutSize = 0;
    unsigned char* tmpByteData;
    if (r2==0)
    {
      SZ_compress_args_int32_NoCkRngeNoGzip_1D(&tmpByteData, oriData, r1, realPrecision, &tmpOutSize, valueRangeSize, minValue);
    }
 
    else
    {
//      printf("Error: doesn't support 5 dimensions for now.\n");
      status = SZ_DERR; //dimension error
    }
    //Call Gzip to do the further compression.
    if(confparams_cpr->szMode==SZ_BEST_SPEED)
    {
      *outSize = tmpOutSize;
      *newByteData = tmpByteData;
    }
    else if(confparams_cpr->szMode==SZ_BEST_COMPRESSION || confparams_cpr->szMode==SZ_DEFAULT_COMPRESSION)
    {
      *outSize = sz_lossless_compress(confparams_cpr->losslessCompressor, confparams_cpr->gzipMode, tmpByteData, tmpOutSize, newByteData);
      free(tmpByteData);
    }
    else
    {
   //   printf("Error: Wrong setting of confparams_cpr->szMode in the int32_t compression.\n");
      status = SZ_MERR; //mode error      
    }
  }
  
  return status;
}

unsigned char* SZ_compress_args(int dataType, void *data, size_t *outSize, int errBoundMode, double absErrBound, 
double relBoundRatio, double pwrBoundRatio, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{

    //TODO
  confparams_cpr->dataType = dataType;
  if(dataType == SZ_FLOAT)
  {
    unsigned char *newByteData = NULL;  
    SZ_compress_args_float(&newByteData, (float *)data, r5, r4, r3, r2, r1, 
    outSize, errBoundMode, absErrBound, relBoundRatio, pwrBoundRatio);  
    return newByteData;
  }
  else if(dataType==SZ_INT8)
  {
    unsigned char *newByteData;
    SZ_compress_args_int8(&newByteData,(int8_t *) data, r5, r4, r3, r2, r1, outSize, errBoundMode, absErrBound, relBoundRatio);
    return newByteData;
  }
  else if(dataType==SZ_INT16)
  {
    unsigned char *newByteData;
    SZ_compress_args_int16(&newByteData,(int16_t *) data, r5, r4, r3, r2, r1, outSize, errBoundMode, absErrBound, relBoundRatio);
    return newByteData;   
  }

  else if(dataType==SZ_INT32) //int type
  {
    unsigned char *newByteData;
    SZ_compress_args_int32(&newByteData, (int32_t *)data, r5, r4, r3, r2, r1, outSize, errBoundMode, absErrBound, relBoundRatio);
    return newByteData;
  }
  else if(dataType==SZ_DOUBLE)
  {
    unsigned char *newByteData;
    SZ_compress_args_double(-1, confparams_cpr->withRegression, &newByteData, (double *)data,r5, r4, r3, r2, r1, 
    outSize, errBoundMode, absErrBound, relBoundRatio, pwrBoundRatio);
    
    return newByteData;
  }
  else
  {
   // printf("Error: dataType can only be SZ_FLOAT, SZ_DOUBLE, SZ_INT8/16/32/64 or SZ_UINT8/16/32/64.\n");
    return NULL;
  }
}


unsigned char *SZ_compress(int dataType, void *data, size_t *outSize, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{ 
  unsigned char *newByteData = SZ_compress_args(dataType, data, outSize, confparams_cpr->errorBoundMode, confparams_cpr->absErrBound, confparams_cpr->relBoundRatio, 
  confparams_cpr->pw_relBoundRatio, r5, r4, r3, r2, r1);
  return newByteData;
}

unsigned long zlib_uncompress5(unsigned char* compressBytes, unsigned long cmpSize, unsigned char** oriData, unsigned long targetOriSize)
{

    int err;
    z_stream d_stream = {0}; /* decompression stream */

    *oriData = (unsigned char*)malloc(sizeof(unsigned char)*targetOriSize);

    d_stream.zalloc = (alloc_func)0;
    d_stream.zfree = (free_func)0;
    d_stream.opaque = (voidpf)0;

    d_stream.next_in  = compressBytes;
    d_stream.avail_in = 0;
    d_stream.next_out = *oriData;
    err = inflateInit(&d_stream);
    while (d_stream.total_out < targetOriSize && d_stream.total_in < cmpSize) {
        d_stream.avail_in = d_stream.avail_out = SZ_ZLIB_BUFFER_SIZE; /* force small buffers */
        err = inflate(&d_stream, Z_SYNC_FLUSH);
        if (err == Z_STREAM_END) break;
    }

    err = inflateEnd(&d_stream);
//    printf("in main szlib uncompress  %d\n",d_stream.total_out);
    return d_stream.total_out;
}

unsigned long sz_lossless_decompress(int losslessCompressor, unsigned char* compressBytes, unsigned long cmpSize, unsigned char** oriData, unsigned long targetOriSize)
{
    unsigned long outSize = 0;
    switch(losslessCompressor)
    {
        case GZIP_COMPRESSOR:
            outSize = zlib_uncompress5(compressBytes, cmpSize, oriData, targetOriSize);
            break;
     
        case ZSTD_COMPRESSOR:
            *oriData = (unsigned char*)malloc(targetOriSize);
            outSize = targetOriSize;
            break;
        default:
           ; //Serial.println("Error: Unrecognized lossless compressor in sz_lossless_decompress()\n");            
    }
    
    return outSize;
}


void new_TightDataPointStorageF_Empty(TightDataPointStorageF **mythis)
{
    *mythis = (TightDataPointStorageF*)malloc(sizeof(TightDataPointStorageF));
    (*mythis)->dataSeriesLength = 0;
    (*mythis)->allSameData = 0;
    (*mythis)->exactDataNum = 0;
    (*mythis)->reservedValue = 0;
    (*mythis)->reqLength = 0;
    (*mythis)->radExpo = 0;

    (*mythis)->rtypeArray = NULL;
    (*mythis)->rtypeArray_size = 0;

    (*mythis)->typeArray = NULL; //its size is dataSeriesLength/4 (or xxx/4+1)
    (*mythis)->typeArray_size = 0;

    (*mythis)->leadNumArray = NULL; //its size is exactDataNum/4 (or exactDataNum/4+1)
    (*mythis)->leadNumArray_size = 0;

    (*mythis)->exactMidBytes = NULL;
    (*mythis)->exactMidBytes_size = 0;

    (*mythis)->residualMidBits = NULL;
    (*mythis)->residualMidBits_size = 0;

    (*mythis)->intervals = 0;
    (*mythis)->isLossless = 0;

    (*mythis)->segment_size = 0;
    (*mythis)->pwrErrBoundBytes = NULL;
    (*mythis)->pwrErrBoundBytes_size = 0;

    (*mythis)->raBytes = NULL;
    (*mythis)->raBytes_size = 0;
}


void new_TightDataPointStorageI_Empty(TightDataPointStorageI **mythis)
{
  *mythis = (TightDataPointStorageI*)malloc(sizeof(TightDataPointStorageI));

  (*mythis)->dataSeriesLength = 0;
  (*mythis)->allSameData = 0;
  (*mythis)->exactDataNum = 0;
  (*mythis)->realPrecision = 0;
  (*mythis)->minValue = 0;
  (*mythis)->exactByteSize = 0;

  (*mythis)->typeArray = NULL; //its size is dataSeriesLength/4 (or xxx/4+1) 
  (*mythis)->typeArray_size = 0;
  
  (*mythis)->exactDataBytes = NULL;
  (*mythis)->exactDataBytes_size = 0;

  (*mythis)->intervals = 0;
  (*mythis)->isLossless = 0;  
}


void new_TightDataPointStorageD_Empty(TightDataPointStorageD **mythis)
{
  *mythis = (TightDataPointStorageD*)malloc(sizeof(TightDataPointStorageD));
  (*mythis)->dataSeriesLength = 0;
  (*mythis)->allSameData = 0;
  (*mythis)->exactDataNum = 0;
  (*mythis)->reservedValue = 0;
  (*mythis)->reqLength = 0;
  (*mythis)->radExpo = 0;

  (*mythis)->rtypeArray = NULL;
  (*mythis)->rtypeArray_size = 0;

  (*mythis)->typeArray = NULL; //its size is dataSeriesLength/4 (or xxx/4+1) 
  (*mythis)->typeArray_size = 0;

  (*mythis)->leadNumArray = NULL; //its size is exactDataNum/4 (or exactDataNum/4+1)
  (*mythis)->leadNumArray_size = 0;

  (*mythis)->exactMidBytes = NULL;
  (*mythis)->exactMidBytes_size = 0;

  (*mythis)->residualMidBits = NULL;
  (*mythis)->residualMidBits_size = 0;
  
  (*mythis)->intervals = 0;
  (*mythis)->isLossless = 0;
  
  (*mythis)->segment_size = 0;
  (*mythis)->pwrErrBoundBytes = NULL;
  (*mythis)->pwrErrBoundBytes_size = 0;
  
  (*mythis)->raBytes = NULL;
  (*mythis)->raBytes_size = 0;

}

short bytesToInt16_bigEndian(unsigned char* bytes)
{
    int temp = 0;
    short res = 0;

    temp = bytes[0] & 0xff;
    res |= temp;

    res <<= 8;
    temp = bytes[1] & 0xff;
    res |= temp;

    return res;
}


float bytesToFloat(unsigned char* bytes)
{
  lfloat buf;
  memcpy(buf.byte, bytes, 4);
  if(sysEndianType==LITTLE_ENDIAN_SYSTEM)
    symTransform_4bytes(buf.byte);
  return buf.value;
}

int bytesToInt(unsigned char* bytes)
{
  lfloat buf;
  memcpy(buf.byte, bytes, 4);
  return buf.ivalue;
}

int bytesToInt32_bigEndian(unsigned char* bytes)
{
  int temp = 0;
  int res = 0;

  res <<= 8;
  temp = bytes[0] & 0xff;
  res |= temp;

  res <<= 8;
  temp = bytes[1] & 0xff;
  res |= temp;

  res <<= 8;
  temp = bytes[2] & 0xff;
  res |= temp;

  res <<= 8;
  temp = bytes[3] & 0xff;
  res |= temp;

  return res;
}



sz_params* convertBytesToSZParams(unsigned char* bytes)
{ 
  sz_params* params = (sz_params*)malloc(sizeof(struct sz_params));
  unsigned char flag1 = bytes[0];
  exe_params->optQuantMode = flag1 >> 7;
  dataEndianType = (flag1 & 0x7f) >> 7;
  sysEndianType = (flag1 & 0x3f) >> 7;

  params->szMode = (flag1 & 0x1f) >> 7;

  int tmp = (flag1 & 0x0f) >> 6;
  switch(tmp)
  {
  case 0:
    params->gzipMode = Z_BEST_SPEED;
    break;
  case 1:
    params->gzipMode = Z_DEFAULT_STRATEGY;
    break;
  case 2:
    params->gzipMode = Z_BEST_COMPRESSION;
    break;
  }

  params->pwr_type = (flag1 & 0x03) >> 6;

    params->sampleDistance = bytesToInt16_bigEndian(&bytes[1]);

  params->predThreshold = 1.0*bytesToInt16_bigEndian(&bytes[3])/10000.0;

    params->dataType = bytes[5] & 0x07;

  params->errorBoundMode = (bytes[5] & 0xf0) >> 4;

    switch(params->errorBoundMode)
    {
  case ABS:
    params->absErrBound = bytesToFloat(&bytes[6]);
    break;
  case REL:
    params->relBoundRatio = bytesToFloat(&bytes[10]);
    break;
  case ABS_AND_REL:
  case ABS_OR_REL:
    params->absErrBound = bytesToFloat(&bytes[6]);
    params->relBoundRatio = bytesToFloat(&bytes[10]);
    break;
  case PSNR:
    params->psnr = bytesToFloat(&bytes[6]);
    break;

  }

    //segment_size  // 2 bytes
     params->segment_size = bytesToInt16_bigEndian(&bytes[14]);

    if(exe_params->optQuantMode==1)
    {
    params->max_quant_intervals = bytesToInt32_bigEndian(&bytes[16]);
    params->quantization_intervals = 0;
  }
  else
  {
    params->max_quant_intervals = 0;
    params->quantization_intervals = bytesToInt32_bigEndian(&bytes[16]);
  }
  return params;
}


int bytesToInt_bigEndian(unsigned char* bytes)
{
    int temp = 0;
    int res = 0;

    res <<= 8;
    temp = bytes[0] & 0xff;
    res |= temp;

    res <<= 8;
    temp = bytes[1] & 0xff;
    res |= temp;

    res <<= 8;
    temp = bytes[2] & 0xff;
    res |= temp;

    res <<= 8;
    temp = bytes[3] & 0xff;
    res |= temp;

    return res;
}


long bytesToLong_bigEndian(unsigned char* b) { 
  long temp = 0;
  long res = 0;

  res <<= 8;
  temp = b[0] & 0xff;
  res |= temp;

  res <<= 8;
  temp = b[1] & 0xff;
  res |= temp;

  res <<= 8;
  temp = b[2] & 0xff;
  res |= temp;

  res <<= 8;
  temp = b[3] & 0xff;
  res |= temp;

  res <<= 8;
  temp = b[4] & 0xff;
  res |= temp;

  res <<= 8;
  temp = b[5] & 0xff;
  res |= temp;

  res <<= 8;
  temp = b[6] & 0xff;
  res |= temp;

  res <<= 8;
  temp = b[7] & 0xff;
  res |= temp;

  return res;
}
 

size_t bytesToSize(unsigned char* bytes)
{
  size_t result = 0;
  if(exe_params->SZ_SIZE_TYPE==4)
    result = bytesToInt_bigEndian(bytes);//4
  else
  result = bytesToLong_bigEndian(bytes);//8
  return result;
}


double bytesToDouble(unsigned char* bytes)
{
  ldouble buf;
  memcpy(buf.byte, bytes, 8);
  if(sysEndianType==LITTLE_ENDIAN_SYSTEM)
    symTransform_8bytes(buf.byte);
  return buf.value;
}


int new_TightDataPointStorageF_fromFlatBytes(TightDataPointStorageF **mythis, unsigned char* flatBytes, size_t flatBytesLength)
{
  new_TightDataPointStorageF_Empty(mythis);
  size_t i, index = 0;
  size_t pwrErrBoundBytes_size = 0, segmentL = 0, radExpoL = 0, pwrErrBoundBytesL = 0;
  char version[3];
  for (i = 0; i < 3; i++)
    version[i] = flatBytes[index++]; //3
  unsigned char sameRByte = flatBytes[index++]; //1
  
  int same = sameRByte & 0x01;
  //confparams_dec->szMode = (sameRByte & 0x06)>>1;
  (*mythis)->isLossless = (sameRByte & 0x10)>>4;
  int isPW_REL = (sameRByte & 0x20)>>5;
  exe_params->SZ_SIZE_TYPE = ((sameRByte & 0x40)>>6)==1?8:4;
  confparams_dec->randomAccess = (sameRByte & 0x02) >> 1;

  int errorBoundMode = ABS;

  sz_params* params = convertBytesToSZParams(&(flatBytes[index]));
  int mode = confparams_dec->szMode;
  int predictionMode = confparams_dec->predictionMode;
  int losslessCompressor = confparams_dec->losslessCompressor;
  int randomAccess = confparams_dec->randomAccess;
  if(confparams_dec!=NULL)
    free(confparams_dec);
  confparams_dec = params;
  confparams_dec->szMode = mode;
  confparams_dec->losslessCompressor = losslessCompressor;
  confparams_dec->randomAccess = randomAccess;

  if(mode==SZ_TEMPORAL_COMPRESSION)
  {
    confparams_dec->szMode = SZ_TEMPORAL_COMPRESSION;
    confparams_dec->predictionMode = predictionMode;
  }

  index += MetaDataByteLength;

  int isRandomAccess = (sameRByte >> 7) & 0x01;

  unsigned char dsLengthBytes[8];
  for (i = 0; i < exe_params->SZ_SIZE_TYPE; i++)
    dsLengthBytes[i] = flatBytes[index++];
  (*mythis)->dataSeriesLength = bytesToSize(dsLengthBytes);// 4 or 8

  if((*mythis)->isLossless==1)
  {
    //(*mymythis)->exactMidBytes = flatBytes+8;
    return errorBoundMode;
  }
  else if(same==1)
  {
    (*mythis)->allSameData = 1;
    size_t exactMidBytesLength = sizeof(float); 
    if(exactMidBytesLength>0)
      (*mythis)->exactMidBytes = (unsigned char*)malloc(sizeof(unsigned char)*exactMidBytesLength);
    else
      (*mythis)->exactMidBytes = NULL;
    for(i = 0;i<exactMidBytesLength;i++)
      (*mythis)->exactMidBytes[i] = flatBytes[index++];
    return errorBoundMode;
  }
  else
    (*mythis)->allSameData = 0;
  if(isRandomAccess == 1)
  {
    (*mythis)->raBytes_size = flatBytesLength - 3 - 1 - MetaDataByteLength - exe_params->SZ_SIZE_TYPE;
    (*mythis)->raBytes = &(flatBytes[index]);
    return errorBoundMode;
  }

  int rtype_ = sameRByte & 0x08;    //=00001000
  unsigned char byteBuf[8];

  for (i = 0; i < 4; i++)
    byteBuf[i] = flatBytes[index++];
  int max_quant_intervals = bytesToInt_bigEndian(byteBuf);// 4

  confparams_dec->maxRangeRadius = max_quant_intervals/2;

  {
    pwrErrBoundBytes_size = 0;
    (*mythis)->pwrErrBoundBytes = NULL;
  }
  for (i = 0; i < 4; i++)
    byteBuf[i] = flatBytes[index++];
  (*mythis)->intervals = bytesToInt_bigEndian(byteBuf);// 4

  for (i = 0; i < 4; i++)
    byteBuf[i] = flatBytes[index++];
  (*mythis)->medianValue = bytesToFloat(byteBuf); //4

  (*mythis)->reqLength = flatBytes[index++]; //1

  for (i = 0; i < 8; i++)
    byteBuf[i] = flatBytes[index++];
  (*mythis)->realPrecision = bytesToDouble(byteBuf);//8

  for (i = 0; i < exe_params->SZ_SIZE_TYPE; i++)
    byteBuf[i] = flatBytes[index++];
  (*mythis)->typeArray_size = bytesToSize(byteBuf);// 4
  if(rtype_!=0)
  {
    for(i = 0;i<exe_params->SZ_SIZE_TYPE;i++)
      byteBuf[i] = flatBytes[index++];
    (*mythis)->rtypeArray_size = bytesToSize(byteBuf);//(ST)
  }
  else
    (*mythis)->rtypeArray_size = 0;

  for (i = 0; i < exe_params->SZ_SIZE_TYPE; i++)
    byteBuf[i] = flatBytes[index++];
  (*mythis)->exactDataNum = bytesToSize(byteBuf);// ST

  for (i = 0; i < exe_params->SZ_SIZE_TYPE; i++)
    byteBuf[i] = flatBytes[index++];
  (*mythis)->exactMidBytes_size = bytesToSize(byteBuf);// ST

  if (rtype_ != 0) {
    if((*mythis)->rtypeArray_size>0)
      (*mythis)->rtypeArray = (unsigned char*)malloc(sizeof(unsigned char)*(*mythis)->rtypeArray_size);
    else
      (*mythis)->rtypeArray = NULL;

    for (i = 0; i < 4; i++)
      byteBuf[i] = flatBytes[index++];
    (*mythis)->reservedValue = bytesToFloat(byteBuf);//4
  }

  size_t logicLeadNumBitsNum = (*mythis)->exactDataNum * 2;
  if (logicLeadNumBitsNum % 8 == 0)
  {
    (*mythis)->leadNumArray_size = logicLeadNumBitsNum >> 3;
  }
  else
  {
    (*mythis)->leadNumArray_size = (logicLeadNumBitsNum >> 3) + 1;
  }

  int minLogValueSize = 0;
  
  if ((*mythis)->rtypeArray != NULL)
  {
    (*mythis)->residualMidBits_size = flatBytesLength - 3 - 1 - MetaDataByteLength - exe_params->SZ_SIZE_TYPE - 4 - radExpoL - segmentL - pwrErrBoundBytesL - 4 - 4 - 1 - 8
        - exe_params->SZ_SIZE_TYPE - exe_params->SZ_SIZE_TYPE - exe_params->SZ_SIZE_TYPE - minLogValueSize - exe_params->SZ_SIZE_TYPE - 4 - (*mythis)->rtypeArray_size
        - minLogValueSize - (*mythis)->typeArray_size - (*mythis)->leadNumArray_size
        - (*mythis)->exactMidBytes_size - pwrErrBoundBytes_size;
    for (i = 0; i < (*mythis)->rtypeArray_size; i++)
      (*mythis)->rtypeArray[i] = flatBytes[index++];
  }
  else
  {
    (*mythis)->residualMidBits_size = flatBytesLength - 3 - 1 - MetaDataByteLength - exe_params->SZ_SIZE_TYPE - 4 - radExpoL - segmentL - pwrErrBoundBytesL - 4 - 4 - 1 - 8
        - exe_params->SZ_SIZE_TYPE - exe_params->SZ_SIZE_TYPE - exe_params->SZ_SIZE_TYPE - minLogValueSize - (*mythis)->typeArray_size
        - (*mythis)->leadNumArray_size - (*mythis)->exactMidBytes_size - pwrErrBoundBytes_size;
  }

  (*mythis)->typeArray = &flatBytes[index];
  //retrieve the number of states (i.e., stateNum)
  (*mythis)->allNodes = bytesToInt_bigEndian((*mythis)->typeArray);    //the first 4 bytes store the stateNum
  (*mythis)->stateNum = ((*mythis)->allNodes+1)/2;

  index+=(*mythis)->typeArray_size;

  (*mythis)->pwrErrBoundBytes = &flatBytes[index];

  index+=pwrErrBoundBytes_size;

  (*mythis)->leadNumArray = &flatBytes[index];

  index+=(*mythis)->leadNumArray_size;

  (*mythis)->exactMidBytes = &flatBytes[index];

  index+=(*mythis)->exactMidBytes_size;

  (*mythis)->residualMidBits = &flatBytes[index];

  //index+=(*mymythis)->residualMidBits_size;

  return errorBoundMode;
}



int new_TightDataPointStorageI_fromFlatBytes(TightDataPointStorageI **mythis, unsigned char* flatBytes, size_t flatBytesLength)
{
  new_TightDataPointStorageI_Empty(mythis);
  size_t i, index = 0;
  char version[3];
  for (i = 0; i < 3; i++)
    version[i] = flatBytes[index++]; //3
  unsigned char sameRByte = flatBytes[index++]; //1

  int same = sameRByte & 0x01;
  //conf_params->szMode = (sameRByte & 0x06)>>1;
  int dataByteSizeCode = (sameRByte & 0x0C)>>2;
  convertDataTypeSizeCode(dataByteSizeCode); //in bytes
  (*mythis)->isLossless = (sameRByte & 0x10)>>4;

  exe_params->SZ_SIZE_TYPE = ((sameRByte & 0x40)>>6)==1?8:4;
  int errorBoundMode = ABS;
  
  if(confparams_dec==NULL)
  {
    confparams_dec = (sz_params*)malloc(sizeof(sz_params));
    memset(confparams_dec, 0, sizeof(sz_params));
  } 
  convertBytesToSZParams(&(flatBytes[index]));
  
  index += MetaDataByteLength; //20 
  
  if(same==0)
    (*mythis)->exactByteSize = flatBytes[index++]; //1
  
  unsigned char dsLengthBytes[8];
  for (i = 0; i < exe_params->SZ_SIZE_TYPE; i++)
    dsLengthBytes[i] = flatBytes[index++];
  (*mythis)->dataSeriesLength = bytesToSize(dsLengthBytes);// ST
  if((*mythis)->isLossless==1)
  {
    //(*mythis)->exactMidBytes = flatBytes+8;
    return errorBoundMode;
  }
  else if(same==1)
  {
    (*mythis)->allSameData = 1;
    (*mythis)->exactDataBytes = &(flatBytes[index]);
    return errorBoundMode;
  }
  else
    (*mythis)->allSameData = 0;

  unsigned char byteBuf[8];

  for (i = 0; i < 4; i++)
    byteBuf[i] = flatBytes[index++];
  int max_quant_intervals = bytesToInt_bigEndian(byteBuf);// 4  

  confparams_dec->maxRangeRadius = max_quant_intervals/2;

  if(errorBoundMode>=PW_REL)
  {
  //  printf("Error: errorBoundMode>=PW_REL in new_TightDataPointStorageI_fromFlatBytes!! Wrong...\n");
    exit(0);
  }

  for (i = 0; i < 4; i++)
    byteBuf[i] = flatBytes[index++];
  (*mythis)->intervals = bytesToInt_bigEndian(byteBuf);// 4 

  for (i = 0; i < 8; i++)
    byteBuf[i] = flatBytes[index++];
  (*mythis)->minValue = bytesToLong_bigEndian(byteBuf); //8
    
  for (i = 0; i < 8; i++)
    byteBuf[i] = flatBytes[index++];
  (*mythis)->realPrecision = bytesToDouble(byteBuf);//8
  
  for (i = 0; i < exe_params->SZ_SIZE_TYPE; i++)
    byteBuf[i] = flatBytes[index++];
  (*mythis)->typeArray_size = bytesToSize(byteBuf);// ST    

  for (i = 0; i < exe_params->SZ_SIZE_TYPE; i++)
    byteBuf[i] = flatBytes[index++];
  (*mythis)->exactDataNum = bytesToSize(byteBuf);// ST
  
  for (i = 0; i < exe_params->SZ_SIZE_TYPE; i++)
    byteBuf[i] = flatBytes[index++];
  (*mythis)->exactDataBytes_size = bytesToSize(byteBuf);// ST   


  (*mythis)->typeArray = &flatBytes[index];
  //retrieve the number of states (i.e., stateNum)
  (*mythis)->allNodes = bytesToInt_bigEndian((*mythis)->typeArray); //the first 4 bytes store the stateNum
  (*mythis)->stateNum = ((*mythis)->allNodes+1)/2;    

  index+=(*mythis)->typeArray_size;
  
  if((*mythis)->exactDataBytes_size > 0)
  { 
    (*mythis)->exactDataBytes = &flatBytes[index];
    index+=(*mythis)->exactDataBytes_size*sizeof(char); 
  }
  else
    (*mythis)->exactDataBytes = NULL; 
  return errorBoundMode;
}


int new_TightDataPointStorageD_fromFlatBytes(TightDataPointStorageD **mythis, unsigned char* flatBytes, size_t flatBytesLength)
{
  new_TightDataPointStorageD_Empty(mythis);
  size_t i, index = 0;
  size_t pwrErrBoundBytes_size = 0, segmentL = 0, radExpoL = 0, pwrErrBoundBytesL = 0;
  char version[3];
  for (i = 0; i < 3; i++)
    version[i] = flatBytes[index++]; //3
  unsigned char sameRByte = flatBytes[index++]; //1
  
  int same = sameRByte & 0x01;
  //confparams_dec->szMode = (sameRByte & 0x06)>>1;
  (*mythis)->isLossless = (sameRByte & 0x10)>>4;
  int isPW_REL = (sameRByte & 0x20)>>5;
  exe_params->SZ_SIZE_TYPE = ((sameRByte & 0x40)>>6)==1?8:4;
  confparams_dec->randomAccess = (sameRByte & 0x02) >> 1;
  
  int errorBoundMode = ABS;  
  if(confparams_dec==NULL)
  {
    confparams_dec = (sz_params*)malloc(sizeof(sz_params));
    memset(confparams_dec, 0, sizeof(sz_params));
  } 
  convertBytesToSZParams(&(flatBytes[index]));

  index += MetaDataByteLength_double;

  int isRegression = (sameRByte >> 7) & 0x01;

  unsigned char dsLengthBytes[8];
  for (i = 0; i < exe_params->SZ_SIZE_TYPE; i++)
    dsLengthBytes[i] = flatBytes[index++];
  (*mythis)->dataSeriesLength = bytesToSize(dsLengthBytes);

  //printf("confparams_dec->szMode=%d\n",confparams_dec->szMode);

  if((*mythis)->isLossless==1)
  {
    //(*this)->exactMidBytes = flatBytes+8;
    return errorBoundMode;
  }
  else if(same==1)
  {
    (*mythis)->allSameData = 1;
    size_t exactMidBytesLength = sizeof(double);//flatBytesLength - 3 - 1 - MetaDataByteLength_double -exe_params->SZ_SIZE_TYPE;
    
    (*mythis)->exactMidBytes = &(flatBytes[index]);
    
    if(exactMidBytesLength>0)
      (*mythis)->exactMidBytes = (unsigned char*)malloc(sizeof(unsigned char)*exactMidBytesLength);
    else
      (*mythis)->exactMidBytes = NULL;
    for(i = 0;i<exactMidBytesLength;i++)
      (*mythis)->exactMidBytes[i] = flatBytes[index++];
    return errorBoundMode;
  
  }
  
  
  else
    (*mythis)->allSameData = 0;
    
  if(isRegression == 1)
  {
    (*mythis)->raBytes_size = flatBytesLength - 3 - 1 - MetaDataByteLength_double - exe_params->SZ_SIZE_TYPE;
    (*mythis)->raBytes = &(flatBytes[index]);
    return errorBoundMode;
  }         
    
  //int rtype_ = 0;//sameRByte & 0x08; //1000   
  int rtype_ = sameRByte & 0x08; //1000   

  unsigned char byteBuf[8];

  for (i = 0; i < 4; i++)
    byteBuf[i] = flatBytes[index++];
  int max_quant_intervals = bytesToInt_bigEndian(byteBuf);// 4  

  confparams_dec->maxRangeRadius = max_quant_intervals/2;

  {
    pwrErrBoundBytes_size = 0;
    (*mythis)->pwrErrBoundBytes = NULL;
  }

  for (i = 0; i < 4; i++)
    byteBuf[i] = flatBytes[index++];
  (*mythis)->intervals = bytesToInt_bigEndian(byteBuf);// 4 

  for (i = 0; i < 8; i++)
    byteBuf[i] = flatBytes[index++];
  (*mythis)->medianValue = bytesToDouble(byteBuf);//8

  (*mythis)->reqLength = flatBytes[index++]; //1
  
  for (i = 0; i < 8; i++)
    byteBuf[i] = flatBytes[index++];
  (*mythis)->realPrecision = bytesToDouble(byteBuf);//8

  for (i = 0; i < exe_params->SZ_SIZE_TYPE; i++)
    byteBuf[i] = flatBytes[index++];
  (*mythis)->typeArray_size = bytesToSize(byteBuf);// exe_params->SZ_SIZE_TYPE  

  if(rtype_!=0)
  {
    for(i = 0;i<exe_params->SZ_SIZE_TYPE;i++) 
      byteBuf[i] = flatBytes[index++];
    (*mythis)->rtypeArray_size = bytesToSize(byteBuf);//ST    
  }
  else
    (*mythis)->rtypeArray_size = 0;

  for (i = 0; i < exe_params->SZ_SIZE_TYPE; i++)
    byteBuf[i] = flatBytes[index++];
  (*mythis)->exactDataNum = bytesToSize(byteBuf);// ST

  for (i = 0; i < exe_params->SZ_SIZE_TYPE; i++)
    byteBuf[i] = flatBytes[index++];
  (*mythis)->exactMidBytes_size = bytesToSize(byteBuf);// ST

  if (rtype_ != 0) {
    if((*mythis)->rtypeArray_size>0)
      (*mythis)->rtypeArray = (unsigned char*)malloc(sizeof(unsigned char)*(*mythis)->rtypeArray_size);
    else
      (*mythis)->rtypeArray = NULL;

    for (i = 0; i < 8; i++)
      byteBuf[i] = flatBytes[index++];
    (*mythis)->reservedValue = bytesToDouble(byteBuf);//8
  }

  size_t logicLeadNumBitsNum = (*mythis)->exactDataNum * 2;
  if (logicLeadNumBitsNum % 8 == 0)
  {
    (*mythis)->leadNumArray_size = logicLeadNumBitsNum >> 3;
  }
  else
  {
    (*mythis)->leadNumArray_size = (logicLeadNumBitsNum >> 3) + 1;
  }
  
  int minLogValueSize = 0;

  if ((*mythis)->rtypeArray != NULL) 
  {
    (*mythis)->residualMidBits_size = flatBytesLength - 3 - 1 - MetaDataByteLength_double - exe_params->SZ_SIZE_TYPE - 4 - radExpoL - segmentL - pwrErrBoundBytesL - 4 - 8 - 1 - 8 
        - exe_params->SZ_SIZE_TYPE - exe_params->SZ_SIZE_TYPE - exe_params->SZ_SIZE_TYPE - minLogValueSize - exe_params->SZ_SIZE_TYPE - 8 - (*mythis)->rtypeArray_size 
        - minLogValueSize - (*mythis)->typeArray_size - (*mythis)->leadNumArray_size
        - (*mythis)->exactMidBytes_size - pwrErrBoundBytes_size - 1 - 1;
    for (i = 0; i < (*mythis)->rtypeArray_size; i++)
      (*mythis)->rtypeArray[i] = flatBytes[index++];
  }
  else
  {
    (*mythis)->residualMidBits_size = flatBytesLength - 3 - 1 - MetaDataByteLength_double - exe_params->SZ_SIZE_TYPE - 4 - radExpoL - segmentL - pwrErrBoundBytesL - 4 - 8 - 1 - 8
        - exe_params->SZ_SIZE_TYPE - exe_params->SZ_SIZE_TYPE - exe_params->SZ_SIZE_TYPE - minLogValueSize - (*mythis)->typeArray_size
        - (*mythis)->leadNumArray_size - (*mythis)->exactMidBytes_size - pwrErrBoundBytes_size - 1 - 1;
  } 

  (*mythis)->typeArray = &flatBytes[index];
  //retrieve the number of states (i.e., stateNum)
  (*mythis)->allNodes = bytesToInt_bigEndian((*mythis)->typeArray); //the first 4 bytes store the stateNum
  (*mythis)->stateNum = ((*mythis)->allNodes+1)/2;  

  index+=(*mythis)->typeArray_size;
  
  (*mythis)->pwrErrBoundBytes = &flatBytes[index];
  
  index+=pwrErrBoundBytes_size;
  
  (*mythis)->leadNumArray = &flatBytes[index];
  
  index+=(*mythis)->leadNumArray_size;
  
  (*mythis)->exactMidBytes = &flatBytes[index];
  
  index+=(*mythis)->exactMidBytes_size;
  
  (*mythis)->residualMidBits = &flatBytes[index];
  
  //index+=(*this)->residualMidBits_size;
  
  return errorBoundMode;
}


int computeDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
    int dimension;
    if(r1==0)
    {
        dimension = 0;
    }
    else if(r2==0)
    {
        dimension = 1;
    }
    else if(r3==0)
    {
        dimension = 2;
    }
    else if(r4==0)
    {
        dimension = 3;
    }
    else if(r5==0)
    {
        dimension = 4;
    }
    else
    {
        dimension = 5;
    }
    return dimension;
}

void convertByteArray2IntArray_fast_2b(size_t stepLength, unsigned char* byteArray, size_t byteArrayLength, unsigned char **intArray)
{
    if(stepLength > byteArrayLength*4)
    {
      //  printf("Error: stepLength > byteArray.length*4\n");
        //printf("stepLength=%zu, byteArray.length=%zu\n", stepLength, byteArrayLength);
        exit(0);
    }
    if(stepLength>0)
        *intArray = (unsigned char*)malloc(stepLength*sizeof(unsigned char));
    else
        *intArray = NULL;
    size_t i, n = 0;

    for (i = 0; i < byteArrayLength; i++) {
        unsigned char tmp = byteArray[i];
        (*intArray)[n++] = (tmp & 0xC0) >> 6;
        if(n==stepLength)
            break;
        (*intArray)[n++] = (tmp & 0x30) >> 4;
        if(n==stepLength)
            break;
        (*intArray)[n++] = (tmp & 0x0C) >> 2;
        if(n==stepLength)
            break;
        (*intArray)[n++] = tmp & 0x03;
        if(n==stepLength)
            break;
    }
}



node new_node2(HuffmanTree *huffmanTree, unsigned int c, unsigned char t)
{
  huffmanTree->pool[huffmanTree->n_nodes].c = c;
  huffmanTree->pool[huffmanTree->n_nodes].t = t;
  return huffmanTree->pool + huffmanTree->n_nodes++;
}
void unpad_tree_uchar(HuffmanTree* huffmanTree, unsigned char* L, unsigned char* R, unsigned int* C, unsigned char *t, unsigned int i, node root)
{
  //root->c = C[i];
  if(root->t==0)
  {
    unsigned char l, r;
    l = L[i];
    if(l!=0)
    {
      node lroot = new_node2(huffmanTree,C[l],t[l]);
      root->left = lroot;
      unpad_tree_uchar(huffmanTree,L,R,C,t,l,lroot);
    }
    r = R[i];
    if(r!=0)
    {
      node rroot = new_node2(huffmanTree,C[r],t[r]);
      root->right = rroot;
      unpad_tree_uchar(huffmanTree,L,R,C,t,r,rroot);
    }
  }
}


void unpad_tree_ushort(HuffmanTree* huffmanTree, unsigned short* L, unsigned short* R, unsigned int* C, unsigned char* t, unsigned int i, node root)
{
  //root->c = C[i];
  if(root->t==0)
  {
    unsigned short l, r;
    l = L[i];
    if(l!=0)
    {
      node lroot = new_node2(huffmanTree,C[l],t[l]);
      root->left = lroot;
      unpad_tree_ushort(huffmanTree,L,R,C,t,l,lroot);
    }
    r = R[i];
    if(r!=0)
    {
      node rroot = new_node2(huffmanTree,C[r],t[r]);
      root->right = rroot;
      unpad_tree_ushort(huffmanTree,L,R,C,t,r,rroot);
    }
  }
}


void unpad_tree_uint(HuffmanTree* huffmanTree, unsigned int* L, unsigned int* R, unsigned int* C, unsigned char* t, unsigned int i, node root)
{
  //root->c = C[i];
  if(root->t==0)
  {
    unsigned int l, r;
    l = L[i];
    if(l!=0)
    {
      node lroot = new_node2(huffmanTree,C[l],t[l]);
      root->left = lroot;
      unpad_tree_uint(huffmanTree,L,R,C,t,l,lroot);
    }
    r = R[i];
    if(r!=0)
    {
      node rroot = new_node2(huffmanTree,C[r],t[r]);
      root->right = rroot;
      unpad_tree_uint(huffmanTree,L,R,C,t,r,rroot);
    }
  }
}



node reconstruct_HuffTree_from_bytes_anyStates(HuffmanTree *huffmanTree, unsigned char* bytes, int nodeCount)
{
  if(nodeCount<=256)
  {
    unsigned char* L = (unsigned char*)malloc(nodeCount*sizeof(unsigned char));
    memset(L, 0, nodeCount*sizeof(unsigned char));
    unsigned char* R = (unsigned char*)malloc(nodeCount*sizeof(unsigned char));
    memset(R, 0, nodeCount*sizeof(unsigned char));
    unsigned int* C = (unsigned int*)malloc(nodeCount*sizeof(unsigned int));
    memset(C, 0, nodeCount*sizeof(unsigned int));
    unsigned char* t = (unsigned char*)malloc(nodeCount*sizeof(unsigned char));
    memset(t, 0, nodeCount*sizeof(unsigned char));
    unsigned char cmpSysEndianType = bytes[0];
    if(cmpSysEndianType!=(unsigned char)sysEndianType)
    {
      unsigned char* p = (unsigned char*)(bytes+1+2*nodeCount*sizeof(unsigned char));
      size_t i = 0, size = nodeCount*sizeof(unsigned int);
      while(1)
      {
        symTransform_4bytes(p);
        i+=sizeof(unsigned int);
        if(i<size)
          p+=sizeof(unsigned int);
        else
          break;
      }
    }
    memcpy(L, bytes+1, nodeCount*sizeof(unsigned char));
    memcpy(R, bytes+1+nodeCount*sizeof(unsigned char), nodeCount*sizeof(unsigned char));
    memcpy(C, bytes+1+2*nodeCount*sizeof(unsigned char), nodeCount*sizeof(unsigned int));
    memcpy(t, bytes+1+2*nodeCount*sizeof(unsigned char)+nodeCount*sizeof(unsigned int), nodeCount*sizeof(unsigned char));
    node root = new_node2(huffmanTree, C[0],t[0]);
    unpad_tree_uchar(huffmanTree,L,R,C,t,0,root);
    free(L);
    free(R);
    free(C);
    free(t);
    return root;
  }
  else if(nodeCount<=65536)
  {
    unsigned short* L = (unsigned short*)malloc(nodeCount*sizeof(unsigned short));
    memset(L, 0, nodeCount*sizeof(unsigned short));
    unsigned short* R = (unsigned short*)malloc(nodeCount*sizeof(unsigned short));
    memset(R, 0, nodeCount*sizeof(unsigned short));
    unsigned int* C = (unsigned int*)malloc(nodeCount*sizeof(unsigned int));
    memset(C, 0, nodeCount*sizeof(unsigned int));
    unsigned char* t = (unsigned char*)malloc(nodeCount*sizeof(unsigned char));
    memset(t, 0, nodeCount*sizeof(unsigned char));

    unsigned char cmpSysEndianType = bytes[0];
    if(cmpSysEndianType!=(unsigned char)sysEndianType)
    {
      unsigned char* p = (unsigned char*)(bytes+1);
      size_t i = 0, size = 3*nodeCount*sizeof(unsigned int);
      while(1)
      {
        symTransform_4bytes(p);
        i+=sizeof(unsigned int);
        if(i<size)
          p+=sizeof(unsigned int);
        else
          break;
      }
    }

    memcpy(L, bytes+1, nodeCount*sizeof(unsigned short));
    memcpy(R, bytes+1+nodeCount*sizeof(unsigned short), nodeCount*sizeof(unsigned short));
    memcpy(C, bytes+1+2*nodeCount*sizeof(unsigned short), nodeCount*sizeof(unsigned int));

    memcpy(t, bytes+1+2*nodeCount*sizeof(unsigned short)+nodeCount*sizeof(unsigned int), nodeCount*sizeof(unsigned char));

    node root = new_node2(huffmanTree,0,0);
    unpad_tree_ushort(huffmanTree,L,R,C,t,0,root);
    free(L);
    free(R);
    free(C);
    free(t);
    return root;
  }
  else //nodeCount>65536
  {
    unsigned int* L = (unsigned int*)malloc(nodeCount*sizeof(unsigned int));
    memset(L, 0, nodeCount*sizeof(unsigned int));
    unsigned int* R = (unsigned int*)malloc(nodeCount*sizeof(unsigned int));
    memset(R, 0, nodeCount*sizeof(unsigned int));
    unsigned int* C = (unsigned int*)malloc(nodeCount*sizeof(unsigned int));
    memset(C, 0, nodeCount*sizeof(unsigned int));
    unsigned char* t = (unsigned char*)malloc(nodeCount*sizeof(unsigned char));
    memset(t, 0, nodeCount*sizeof(unsigned char));
    unsigned char cmpSysEndianType = bytes[0];
    if(cmpSysEndianType!=(unsigned char)sysEndianType)
    {
      unsigned char* p = (unsigned char*)(bytes+1);
      size_t i = 0, size = 3*nodeCount*sizeof(unsigned int);
      while(1)
      {
        symTransform_4bytes(p);
        i+=sizeof(unsigned int);
        if(i<size)
          p+=sizeof(unsigned int);
        else
          break;
      }
    }

    memcpy(L, bytes+1, nodeCount*sizeof(unsigned int));
    memcpy(R, bytes+1+nodeCount*sizeof(unsigned int), nodeCount*sizeof(unsigned int));
    memcpy(C, bytes+1+2*nodeCount*sizeof(unsigned int), nodeCount*sizeof(unsigned int));

    memcpy(t, bytes+1+3*nodeCount*sizeof(unsigned int), nodeCount*sizeof(unsigned char));

    node root = new_node2(huffmanTree,0,0);
    unpad_tree_uint(huffmanTree,L,R,C,t,0,root);
    free(L);
    free(R);
    free(C);
    free(t);
    return root;
  }
}
void decode(unsigned char *s, size_t targetLength, node t, int *out)
{
  size_t i = 0, byteIndex = 0, count = 0;
  int r;
  node n = t;

  if(n->t) //root->t==1 means that all state values are the same (constant)
  {
    for(count=0;count<targetLength;count++)
      out[count] = n->c;
    return;
  }

  for(i=0;count<targetLength;i++)
  {

    byteIndex = i>>3; //i/8
    r = i%8;
    if(((s[byteIndex] >> (7-r)) & 0x01) == 0)
      n = n->left;
    else
      n = n->right;

    if (n->t) {
      //putchar(n->c);
      out[count] = n->c;
      n = t;
      count++;
    }
  }
//  putchar('\n');
  if (t != n) //printf("garbage input\n");
  return;
}

void decode_withTree(HuffmanTree* huffmanTree, unsigned char *s, size_t targetLength, int *out)
{
  size_t encodeStartIndex;
  size_t nodeCount = bytesToInt_bigEndian(s);
  node root = reconstruct_HuffTree_from_bytes_anyStates(huffmanTree,s+8, nodeCount);

  if(nodeCount<=256)
    encodeStartIndex = 1+3*nodeCount*sizeof(unsigned char)+nodeCount*sizeof(unsigned int);
  else if(nodeCount<=65536)
    encodeStartIndex = 1+2*nodeCount*sizeof(unsigned short)+nodeCount*sizeof(unsigned char)+nodeCount*sizeof(unsigned int);
  else
    encodeStartIndex = 1+3*nodeCount*sizeof(unsigned int)+nodeCount*sizeof(unsigned char);
  decode(s+8+encodeStartIndex, targetLength, root, out);
}
 int getRightMovingSteps(int kMod8, int resiBitLength) {
  return 8 - kMod8 - resiBitLength;
}

int getMaskRightCode(int m) {

  switch (m) {
  case 1:
    return 0x01;
  case 2:
    return 0x03;
  case 3:
    return 0x07;
  case 4:
    return 0x0F;
  case 5:
    return 0x1F;
  case 6:
    return 0x3F;
  case 7:
    return 0X7F;
  case 8:
    return 0XFF;
  default:
    return 0;
  }
}


int getRightMovingCode(int kMod8, int resiBitLength)
{
  int rightMovingSteps = 8 - kMod8 - resiBitLength;
  if(rightMovingSteps < 0)
  {
    switch(-rightMovingSteps)
    {
    case 1:
      return 0x80;
    case 2:
      return 0xC0;
    case 3:
      return 0xE0;
    case 4:
      return 0xF0;
    case 5:
      return 0xF8;
    case 6:
      return 0xFC;
    case 7:
      return 0XFE;
    default:
      return 0;
    }
  }
  else //if(rightMovingSteps >= 0)
  {
    int a = getMaskRightCode(8 - kMod8);
    int b = getMaskRightCode(8 - kMod8 - resiBitLength);
    int c = a - b;
    return c;
  }
}
int getLeftMovingCode(int kMod8)
{
  return getMaskRightCode(8 - kMod8);
}


void decompressDataSeries_float_1D(float** data, size_t dataSeriesLength, TightDataPointStorageF* tdps)
{
  updateQuantizationInfo(tdps->intervals);
  size_t i, j, k = 0, p = 0, l = 0; // k is to track the location of residual_bit
                // in resiMidBits, p is to track the
                // byte_index of resiMidBits, l is for
                // leadNum
  unsigned char* leadNum;
  double interval = tdps->realPrecision*2;
  
  convertByteArray2IntArray_fast_2b(tdps->exactDataNum, tdps->leadNumArray, tdps->leadNumArray_size, &leadNum);

  *data = (float*)malloc(sizeof(float)*dataSeriesLength);

  int* type = (int*)malloc(dataSeriesLength*sizeof(int));

  HuffmanTree* huffmanTree = createHuffmanTree(tdps->stateNum);
  decode_withTree(huffmanTree, tdps->typeArray, dataSeriesLength, type);
  SZ_ReleaseHuffman(huffmanTree);

  unsigned char preBytes[4];
  unsigned char curBytes[4];

  memset(preBytes, 0, 4);

  size_t curByteIndex = 0;
  int reqBytesLength, resiBitsLength, resiBits;
  unsigned char leadingNum;
  float medianValue, exactData, predValue;

  reqBytesLength = tdps->reqLength/8;
  resiBitsLength = tdps->reqLength%8;
  //printf("resi %d\n",tdps->reqLength);
  medianValue = tdps->medianValue;

  int type_;
  for (i = 0; i < dataSeriesLength; i++) {
    type_ = type[i];
    switch (type_) {
    case 0:
      // compute resiBits
      resiBits = 0;
      if (resiBitsLength != 0) {
        int kMod8 = k % 8;
        int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
        if (rightMovSteps > 0) {
      int code = getRightMovingCode(kMod8, resiBitsLength);
      resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
        } else if (rightMovSteps < 0) {
          int code1 = getLeftMovingCode(kMod8);
          int code2 = getRightMovingCode(kMod8, resiBitsLength);
          int leftMovSteps = -rightMovSteps;
          rightMovSteps = 8 - leftMovSteps;
          resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
          p++;
          resiBits = resiBits
              | ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
        } else // rightMovSteps == 0

        {
          int code = getRightMovingCode(kMod8, resiBitsLength);
          resiBits = (tdps->residualMidBits[p] & code);
          p++;
        }
             
        k += resiBitsLength;
      }

      // recover the exact data
      memset(curBytes, 0, 4);
      leadingNum = leadNum[l++];
      memcpy(curBytes, preBytes, leadingNum);
      for (j = leadingNum; j < reqBytesLength; j++)
        curBytes[j] = tdps->exactMidBytes[curByteIndex++];
      if (resiBitsLength != 0) {
        unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
        curBytes[reqBytesLength] = resiByte;
      }

      exactData = bytesToFloat(curBytes);
      (*data)[i] = exactData + medianValue;
      memcpy(preBytes,curBytes,4);
      //printf("XXXXXXXXXXXXXX\n");
      break;
    default:
      //predValue = 2 * (*data)[i-1] - (*data)[i-2];
      predValue = (*data)[i-1];
      (*data)[i] = predValue + (type_-exe_params->intvRadius)*interval;
      break;
   }
    //printf("BINGO %.30G\n",(*data)[i]);
  }
  free(leadNum);
  free(type);
  return;
}


void decompressDataSeries_double_1D(double** data, size_t dataSeriesLength, double* hist_data, TightDataPointStorageD* tdps) 
{
  updateQuantizationInfo(tdps->intervals);
  int intvRadius = tdps->intervals/2;
  size_t i, j, k = 0, p = 0, l = 0; // k is to track the location of residual_bit
                // in resiMidBits, p is to track the
                // byte_index of resiMidBits, l is for
                // leadNum
  unsigned char* leadNum;
  double interval = tdps->realPrecision*2;
  
  convertByteArray2IntArray_fast_2b(tdps->exactDataNum, tdps->leadNumArray, tdps->leadNumArray_size, &leadNum);
  *data = (double*)malloc(sizeof(double)*dataSeriesLength);

  int* type = (int*)malloc(dataSeriesLength*sizeof(int));

  HuffmanTree* huffmanTree = createHuffmanTree(tdps->stateNum);
  decode_withTree(huffmanTree, tdps->typeArray, dataSeriesLength, type);
  SZ_ReleaseHuffman(huffmanTree); 
  
  unsigned char preBytes[8];
  unsigned char curBytes[8];
  
  memset(preBytes, 0, 8);

  size_t curByteIndex = 0;
  int reqBytesLength, resiBitsLength, resiBits; 
  unsigned char leadingNum; 
  double medianValue, exactData, predValue;
  
  reqBytesLength = tdps->reqLength/8;
  resiBitsLength = tdps->reqLength%8;
  medianValue = tdps->medianValue;
  
  int type_;
  for (i = 0; i < dataSeriesLength; i++) {
    type_ = type[i];
    switch (type_) {
    case 0:
      // compute resiBits
      resiBits = 0;
      if (resiBitsLength != 0) {
        int kMod8 = k % 8;
        int rightMovSteps = getRightMovingSteps(kMod8, resiBitsLength);
        if (rightMovSteps > 0) {
          int code = getRightMovingCode(kMod8, resiBitsLength);
          resiBits = (tdps->residualMidBits[p] & code) >> rightMovSteps;
        } else if (rightMovSteps < 0) {
          int code1 = getLeftMovingCode(kMod8);
          int code2 = getRightMovingCode(kMod8, resiBitsLength);
          int leftMovSteps = -rightMovSteps;
          rightMovSteps = 8 - leftMovSteps;
          resiBits = (tdps->residualMidBits[p] & code1) << leftMovSteps;
          p++;
          resiBits = resiBits
              | ((tdps->residualMidBits[p] & code2) >> rightMovSteps);
        } else // rightMovSteps == 0
        {
          int code = getRightMovingCode(kMod8, resiBitsLength);
          resiBits = (tdps->residualMidBits[p] & code);
          p++;
        }
        k += resiBitsLength;
      }

      // recover the exact data
      memset(curBytes, 0, 8);
      leadingNum = leadNum[l++];
      memcpy(curBytes, preBytes, leadingNum);
      for (j = leadingNum; j < reqBytesLength; j++)
        curBytes[j] = tdps->exactMidBytes[curByteIndex++];
      if (resiBitsLength != 0) {
        unsigned char resiByte = (unsigned char) (resiBits << (8 - resiBitsLength));
        curBytes[reqBytesLength] = resiByte;
      }
      
      exactData = bytesToDouble(curBytes);
      (*data)[i] = exactData + medianValue;
      memcpy(preBytes,curBytes,8);
      break;
    default:
      //predValue = 2 * (*data)[i-1] - (*data)[i-2];
      predValue = (*data)[i-1];
      (*data)[i] = predValue + (type_-exe_params->intvRadius)*interval;
      break;
    }
  //  printf("%.30G\n",(*data)[i]);
  }

  free(leadNum);
  free(type);
  return;
}


int computeRightShiftBits(int exactByteSize, int dataType)
{
  int rightShift = 0; 
  switch(dataType)
  {
  case SZ_INT8:
    rightShift = 8 - exactByteSize*8;
    break;
    case SZ_INT16:
    rightShift = 16 - exactByteSize*8;
    break;
  case SZ_INT32:
    rightShift = 32 - exactByteSize*8;
    break;
  }
  return rightShift;  
}

void decompressDataSeries_int8_1D(int8_t** data, size_t dataSeriesLength, TightDataPointStorageI* tdps) 
{
  updateQuantizationInfo(tdps->intervals);
  double interval = tdps->realPrecision*2;
  
  *data = (int8_t*)malloc(sizeof(int8_t)*dataSeriesLength);

  int* type = (int*)malloc(dataSeriesLength*sizeof(int));

  HuffmanTree* huffmanTree = createHuffmanTree(tdps->stateNum);
  decode_withTree(huffmanTree, tdps->typeArray, dataSeriesLength, type);
  SZ_ReleaseHuffman(huffmanTree); 

  long predValue, tmp;
  int8_t minValue, exactData;
  
  minValue = tdps->minValue;
  
  int exactByteSize = tdps->exactByteSize;
  unsigned char* exactDataBytePointer = tdps->exactDataBytes;
  
  unsigned char curBytes[8] = {0,0,0,0,0,0,0,0};
  
  int rightShiftBits = computeRightShiftBits(exactByteSize, SZ_INT8);
  if(rightShiftBits<0)
  {
//    printf("Error: rightShift < 0!\n");
    exit(0);
  }
  int type_;
  for (size_t i = 0; i < dataSeriesLength; i++) {
    type_ = type[i];
    switch (type_) {
    case 0:
      // recover the exact data 
      memcpy(curBytes, exactDataBytePointer, exactByteSize);
      exactData = curBytes[0];
      exactData = (uint8_t)exactData >> rightShiftBits;
      exactDataBytePointer += exactByteSize;
      (*data)[i] = exactData + minValue;
      break;
    default:
      //predValue = 2 * (*data)[i-1] - (*data)[i-2];
      predValue = (*data)[i-1];
      tmp = predValue + (type_-exe_params->intvRadius)*interval;
      if(tmp >= SZ_INT8_MIN&&tmp<SZ_INT8_MAX)
        (*data)[i] = tmp;
      else if(tmp < SZ_INT8_MIN)
        (*data)[i] = SZ_INT8_MIN;
      else
        (*data)[i] = SZ_INT8_MAX;
      break;
    }
    //printf("%.30G\n",(*data)[i]);
  }
  free(type);
  return;
}

void decompressDataSeries_int16_1D(int16_t** data, size_t dataSeriesLength, TightDataPointStorageI* tdps) 
{
  updateQuantizationInfo(tdps->intervals);
  size_t i;
  double interval = tdps->realPrecision*2;
  
  *data = (int16_t*)malloc(sizeof(int16_t)*dataSeriesLength);

  int* type = (int*)malloc(dataSeriesLength*sizeof(int));

  HuffmanTree* huffmanTree = createHuffmanTree(tdps->stateNum);
  decode_withTree(huffmanTree, tdps->typeArray, dataSeriesLength, type);
  SZ_ReleaseHuffman(huffmanTree); 
  long predValue, tmp;
  int16_t minValue, exactData;
  
  minValue = tdps->minValue;
  
  int exactByteSize = tdps->exactByteSize;
  unsigned char* exactDataBytePointer = tdps->exactDataBytes;
  
  unsigned char curBytes[8] = {0,0,0,0,0,0,0,0};
  
  int rightShiftBits = computeRightShiftBits(exactByteSize, SZ_INT16);
  if(rightShiftBits<0)
  {
   // printf("Error: rightShift < 0!\n");
    exit(0);
  }
  int type_;
  for (i = 0; i < dataSeriesLength; i++) {
    type_ = type[i];
    switch (type_) {
    case 0:
      // recover the exact data 
      memcpy(curBytes, exactDataBytePointer, exactByteSize);
      exactData = bytesToInt16_bigEndian(curBytes);
      exactData = (uint16_t)exactData >> rightShiftBits;
      exactDataBytePointer += exactByteSize;
      (*data)[i] = exactData + minValue;
      break;
    default:
      //predValue = 2 * (*data)[i-1] - (*data)[i-2];
      predValue = (*data)[i-1];
      tmp = predValue + (type_-exe_params->intvRadius)*interval;
      if(tmp >= SZ_INT16_MIN&&tmp<SZ_INT16_MAX)
        (*data)[i] = tmp;
      else if(tmp < SZ_INT16_MIN)
        (*data)[i] = SZ_INT16_MIN;
      else
        (*data)[i] = SZ_INT16_MAX;
      break;
    }
    //printf("%.30G\n",(*data)[i]);
  }
  free(type);
  return;
}


void decompressDataSeries_int32_1D(int32_t** data, size_t dataSeriesLength, TightDataPointStorageI* tdps) 
{
  updateQuantizationInfo(tdps->intervals);
  size_t i;
  double interval = tdps->realPrecision*2;
  
  *data = (int32_t*)malloc(sizeof(int32_t)*dataSeriesLength);

  int* type = (int*)malloc(dataSeriesLength*sizeof(int));

  HuffmanTree* huffmanTree = createHuffmanTree(tdps->stateNum);
  decode_withTree(huffmanTree, tdps->typeArray, dataSeriesLength, type);
  SZ_ReleaseHuffman(huffmanTree); 
  int32_t minValue, exactData, predValue;
  
  minValue = tdps->minValue;
  
  int exactByteSize = tdps->exactByteSize;
  unsigned char* exactDataBytePointer = tdps->exactDataBytes;
  
  unsigned char curBytes[8] = {0,0,0,0,0,0,0,0};
  
  int rightShiftBits = computeRightShiftBits(exactByteSize, SZ_INT32);
  if(rightShiftBits<0)
  {
  //  printf("Error: rightShift < 0!\n");
    exit(0);
  }
  int type_;
  for (i = 0; i < dataSeriesLength; i++) {
    type_ = type[i];
    switch (type_) {
    case 0:
      // recover the exact data 
      memcpy(curBytes, exactDataBytePointer, exactByteSize);
      exactData = bytesToInt32_bigEndian(curBytes);
      exactData = (uint32_t)exactData >> rightShiftBits;
      exactDataBytePointer += exactByteSize;
      (*data)[i] = exactData + minValue;
      break;
    default:
      //predValue = 2 * (*data)[i-1] - (*data)[i-2];
      predValue = (*data)[i-1];
      (*data)[i] = predValue + (type_-exe_params->intvRadius)*interval;
      break;
    }
    //printf("%.30G\n",(*data)[i]);
  }
  free(type);
  return;
}


char numberOfLeadingZeros_Int(int i) {
  if (i == 0)
    return 32;
  unsigned char n = 1;
  if (((unsigned int)i) >> 16 == 0) { n += 16; i <<= 16; }
  if (((unsigned int)i) >> 24 == 0) { n +=  8; i <<=  8; }
  if (((unsigned int)i) >> 28 == 0) { n +=  4; i <<=  4; }
  if (((unsigned int)i) >> 30 == 0) { n +=  2; i <<=  2; }
  n -= ((unsigned int)i) >> 31;
  return n;
}
char numberOfLeadingZeros_Long(long i) {
   if (i == 0)
    return 64;
  unsigned char n = 1;
  int x = (int)(((unsigned long)i) >> 32);
  if (x == 0) { n += 32; x = (int)i; }
  if (((unsigned int)x) >> 16 == 0) { n += 16; x <<= 16; }
  if (((unsigned int)x) >> 24 == 0) { n +=  8; x <<=  8; }
  if (((unsigned int)x) >> 28 == 0) { n +=  4; x <<=  4; }
  if (((unsigned int)x) >> 30 == 0) { n +=  2; x <<=  2; }
  n -= ((unsigned int)x) >> 31;
  return n;
}

int computeBitNumRequired(size_t dataLength)
{
  if(exe_params->SZ_SIZE_TYPE==4)
    return 32 - numberOfLeadingZeros_Int(dataLength);
  else
    return 64 - numberOfLeadingZeros_Long(dataLength);

}

void getSnapshotData_float_1D(float** data, size_t dataSeriesLength, TightDataPointStorageF* tdps, int errBoundMode)
{
    size_t i;
    if (tdps->allSameData) {
        float value = bytesToFloat(tdps->exactMidBytes);
        *data = (float*)malloc(sizeof(float)*dataSeriesLength);
        for (i = 0; i < dataSeriesLength; i++)
            (*data)[i] = value;
    } else {
        if (tdps->rtypeArray == NULL) {
            if(errBoundMode < PW_REL)
            {
            decompressDataSeries_float_1D(data, dataSeriesLength, tdps);
            }
          
            return;
        } 
    }
}




void getSnapshotData_int8_1D(int8_t** data, size_t dataSeriesLength, TightDataPointStorageI* tdps, int errBoundMode)
{  
  size_t i;

  if (tdps->allSameData) {
    int8_t value = tdps->exactDataBytes[0];
    *data = (int8_t*)malloc(sizeof(int8_t)*dataSeriesLength);
    for (i = 0; i < dataSeriesLength; i++)
      (*data)[i] = value;
  } else {
    decompressDataSeries_int8_1D(data, dataSeriesLength, tdps);
  }
}

void getSnapshotData_int16_1D(int16_t** data, size_t dataSeriesLength, TightDataPointStorageI* tdps, int errBoundMode)
{  
  size_t i;

  if (tdps->allSameData) {
    int16_t value = bytesToInt16_bigEndian(tdps->exactDataBytes);
    *data = (int16_t*)malloc(sizeof(int16_t)*dataSeriesLength);
    for (i = 0; i < dataSeriesLength; i++)
      (*data)[i] = value;
  } else {
    decompressDataSeries_int16_1D(data, dataSeriesLength, tdps);
  }
}


void getSnapshotData_double_1D(double** data, size_t dataSeriesLength, TightDataPointStorageD* tdps, int errBoundMode, int compressionType, double* hist_data) 
{
  size_t i;
  if (tdps->allSameData) {
    double value = bytesToDouble(tdps->exactMidBytes);
    *data = (double*)malloc(sizeof(double)*dataSeriesLength);
    for (i = 0; i < dataSeriesLength; i++)
      (*data)[i] = value;
     
  } else {
    if (tdps->rtypeArray == NULL) {
      if(errBoundMode < PW_REL)
      {
         decompressDataSeries_double_1D(data, dataSeriesLength, hist_data, tdps);
      }
      
      return;
    } 
  }
}


void free_TightDataPointStorageF2(TightDataPointStorageF *tdps)
{
  free(tdps);
}

void free_TightDataPointStorageD2(TightDataPointStorageD *tdps)
{      
  free(tdps);
}
 
int SZ_decompress_args_float(float** newData, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, unsigned char* cmpBytes, size_t cmpSize)
{
  int status = SZ_SCES;
  size_t dataLength = computeDataLength(r5,r4,r3,r2,r1);

  //unsigned char* tmpBytes;
  size_t targetUncompressSize = dataLength <<2; //i.e., *4
  //tmpSize must be "much" smaller than dataLength
  size_t i, tmpSize = 8+MetaDataByteLength+exe_params->SZ_SIZE_TYPE;
  unsigned char* szTmpBytes;

  if(cmpSize!=8+4+MetaDataByteLength && cmpSize!=8+8+MetaDataByteLength) //4,8 means two posibilities of SZ_SIZE_TYPE
  {
    confparams_dec->losslessCompressor = 0;//is_lossless_compressed_data(cmpBytes, cmpSize);
    if(confparams_dec->szMode!=SZ_TEMPORAL_COMPRESSION)
    {
      if(confparams_dec->losslessCompressor!=-1)
        confparams_dec->szMode = SZ_BEST_COMPRESSION;
      else
        confparams_dec->szMode = SZ_BEST_SPEED;
    }

    if(confparams_dec->szMode==SZ_BEST_SPEED)
    {
      tmpSize = cmpSize;
      szTmpBytes = cmpBytes;
    }
    else if(confparams_dec->szMode==SZ_BEST_COMPRESSION || confparams_dec->szMode==SZ_DEFAULT_COMPRESSION || confparams_dec->szMode==SZ_TEMPORAL_COMPRESSION)
    {
      if(targetUncompressSize<MIN_ZLIB_DEC_ALLOMEM_BYTES) //Considering the minimum size
        targetUncompressSize = MIN_ZLIB_DEC_ALLOMEM_BYTES;
      tmpSize = sz_lossless_decompress(confparams_dec->losslessCompressor, cmpBytes, (unsigned long)cmpSize, &szTmpBytes, (unsigned long)targetUncompressSize+4+MetaDataByteLength+exe_params->SZ_SIZE_TYPE);//   (unsigned long)targetUncompressSize+8: consider the total length under lossless compression mode is actually 3+4+1+targetUncompressSize
      
    }
    else
    {
      //printf("Wrong value of confparams_dec->szMode in the double compressed bytes.\n");
      status = SZ_MERR;
      return status;
    }
  }
  else
    szTmpBytes = cmpBytes;
  //TODO: convert szTmpBytes to data array.
  TightDataPointStorageF* tdps;
  int errBoundMode = new_TightDataPointStorageF_fromFlatBytes(&tdps, szTmpBytes, tmpSize);
  int dim = computeDimension(r5,r4,r3,r2,r1);
  int floatSize = sizeof(float);
  if(tdps->isLossless)
  {
    *newData = (float*)malloc(floatSize*dataLength);
    if(sysEndianType==BIG_ENDIAN_SYSTEM)
    {
      memcpy(*newData, szTmpBytes+4+MetaDataByteLength+exe_params->SZ_SIZE_TYPE, dataLength*floatSize);
    }
    else
    {
      unsigned char* p = szTmpBytes+4+MetaDataByteLength+exe_params->SZ_SIZE_TYPE;
      for(i=0;i<dataLength;i++,p+=floatSize)
        (*newData)[i] = bytesToFloat(p);
    }
  }
  else
  {

    if(tdps->raBytes_size > 0) //v2.0
    {            
      if (dim == 1)
        getSnapshotData_float_1D(newData,r1,tdps, errBoundMode);
      
      else
      {
     //   printf("Error: currently support only at most 4 dimensions!\n");
        status = SZ_DERR;
      }
    }
    else //1.4.13
    {
   if (dim == 1)
        getSnapshotData_float_1D(newData,r1,tdps, errBoundMode);
         
      else
      {
     //   printf("Error: currently support only at most 4 dimensions!\n");
        status = SZ_DERR;
      }
    }
  }
  free_TightDataPointStorageF2(tdps);
  if(confparams_dec->szMode!=SZ_BEST_SPEED && cmpSize!=8+MetaDataByteLength+exe_params->SZ_SIZE_TYPE)
    free(szTmpBytes);

  return status;
}


int SZ_decompress_args_int8(int8_t** newData, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, unsigned char* cmpBytes, size_t cmpSize)
{
  int status = SZ_SCES;
  size_t dataLength = computeDataLength(r5,r4,r3,r2,r1);
  
  //unsigned char* tmpBytes;
  size_t targetUncompressSize = dataLength <<2; //i.e., *4
  //tmpSize must be "much" smaller than dataLength
  size_t i, tmpSize = 3+MetaDataByteLength+1+sizeof(int8_t)+exe_params->SZ_SIZE_TYPE;
  unsigned char* szTmpBytes;  
    
  if(cmpSize!=4+1+4+MetaDataByteLength && cmpSize!=4+1+8+MetaDataByteLength)
  {
    confparams_dec->losslessCompressor = 0; //is_lossless_compressed_data(cmpBytes, cmpSize);
    if(confparams_dec->losslessCompressor!=-1)
      confparams_dec->szMode = SZ_BEST_COMPRESSION;
    else
      confparams_dec->szMode = SZ_BEST_SPEED;   
    if(confparams_dec->szMode==SZ_BEST_SPEED)
    {
      tmpSize = cmpSize;
      szTmpBytes = cmpBytes;  
    }
    else if(confparams_dec->szMode==SZ_BEST_COMPRESSION || confparams_dec->szMode==SZ_DEFAULT_COMPRESSION)
    {
      if(targetUncompressSize<MIN_ZLIB_DEC_ALLOMEM_BYTES) //Considering the minimum size
        targetUncompressSize = MIN_ZLIB_DEC_ALLOMEM_BYTES; 
      tmpSize = sz_lossless_decompress(confparams_dec->losslessCompressor, cmpBytes, (unsigned long)cmpSize, &szTmpBytes, (unsigned long)targetUncompressSize+4+MetaDataByteLength+exe_params->SZ_SIZE_TYPE);//   (unsigned long)targetUncompressSize+8: consider the total length under lossless compression mode is actually 3+4+1+targetUncompressSize
      //szTmpBytes = (unsigned char*)malloc(sizeof(unsigned char)*tmpSize);
      //memcpy(szTmpBytes, tmpBytes, tmpSize);
      //free(tmpBytes); //release useless memory    
    }
    else
    {
    //  printf("Wrong value of confparams_dec->szMode in the double compressed bytes.\n");
      status = SZ_MERR;
      return status;
    } 
  }
  else
    szTmpBytes = cmpBytes;
  //TODO: convert szTmpBytes to data array.
  TightDataPointStorageI* tdps;
  int errBoundMode = new_TightDataPointStorageI_fromFlatBytes(&tdps, szTmpBytes, tmpSize);
  //writeByteData(tdps->typeArray, tdps->typeArray_size, "decompress-typebytes.tbt");
  int dim = computeDimension(r5,r4,r3,r2,r1); 
  int intSize = sizeof(int8_t);
  if(tdps->isLossless)
  {
    *newData = (int8_t*)malloc(intSize*dataLength);
    if(sysEndianType==BIG_ENDIAN_SYSTEM)
    {
      memcpy(*newData, szTmpBytes+4+MetaDataByteLength+exe_params->SZ_SIZE_TYPE, dataLength*intSize);
    }
    else
    {
      unsigned char* p = szTmpBytes+4+MetaDataByteLength+exe_params->SZ_SIZE_TYPE;
      for(i=0;i<dataLength;i++,p+=intSize)
        (*newData)[i] = *p;
    }   
  }
  else if (dim == 1)
    getSnapshotData_int8_1D(newData,r1,tdps, errBoundMode);
  else
  {
  //  printf("Error: currently support only at most 4 dimensions!\n");
    status = SZ_DERR;
  }
  free_TightDataPointStorageI2(tdps);
  if(confparams_dec->szMode!=SZ_BEST_SPEED && cmpSize!=4+sizeof(int8_t)+exe_params->SZ_SIZE_TYPE+MetaDataByteLength)
    free(szTmpBytes);
  return status;
}

int SZ_decompress_args_int16(int16_t** newData, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, unsigned char* cmpBytes, size_t cmpSize)
{
  int status = SZ_SCES;
  size_t dataLength = computeDataLength(r5,r4,r3,r2,r1);
  
  //unsigned char* tmpBytes;
  size_t targetUncompressSize = dataLength <<2; //i.e., *4
  //tmpSize must be "much" smaller than dataLength
  size_t i, tmpSize = 3+MetaDataByteLength+1+sizeof(int16_t)+exe_params->SZ_SIZE_TYPE;
  unsigned char* szTmpBytes;  
    
  if(cmpSize!=4+2+4+MetaDataByteLength && cmpSize!=4+2+8+MetaDataByteLength)
  {
    confparams_dec->losslessCompressor = 0 ;  //is_lossless_compressed_data(cmpBytes, cmpSize);
    if(confparams_dec->losslessCompressor!=-1)
      confparams_dec->szMode = SZ_BEST_COMPRESSION;
    else
      confparams_dec->szMode = SZ_BEST_SPEED;   
    if(confparams_dec->szMode==SZ_BEST_SPEED)
    {
      tmpSize = cmpSize;
      szTmpBytes = cmpBytes;  
    }
    else if(confparams_dec->szMode==SZ_BEST_COMPRESSION || confparams_dec->szMode==SZ_DEFAULT_COMPRESSION)
    {
      if(targetUncompressSize<MIN_ZLIB_DEC_ALLOMEM_BYTES) //Considering the minimum size
        targetUncompressSize = MIN_ZLIB_DEC_ALLOMEM_BYTES; 
      tmpSize = sz_lossless_decompress(confparams_dec->losslessCompressor, cmpBytes, (unsigned long)cmpSize, &szTmpBytes, (unsigned long)targetUncompressSize+4+MetaDataByteLength+exe_params->SZ_SIZE_TYPE);//   (unsigned long)targetUncompressSize+8: consider the total length under lossless compression mode is actually 3+4+1+targetUncompressSize
      //szTmpBytes = (unsigned char*)malloc(sizeof(unsigned char)*tmpSize);
      //memcpy(szTmpBytes, tmpBytes, tmpSize);
      //free(tmpBytes); //release useless memory    
    }
    else
    {
     // printf("Wrong value of confparams_dec->szMode in the double compressed bytes.\n");
      status = SZ_MERR;
      return status;
    } 
  }
  else
    szTmpBytes = cmpBytes;
  //TODO: convert szTmpBytes to data array.
  TightDataPointStorageI* tdps;
  int errBoundMode = new_TightDataPointStorageI_fromFlatBytes(&tdps, szTmpBytes, tmpSize);
  //writeByteData(tdps->typeArray, tdps->typeArray_size, "decompress-typebytes.tbt");
  int dim = computeDimension(r5,r4,r3,r2,r1); 
  int intSize = sizeof(int16_t);
  if(tdps->isLossless)
  {
    *newData = (int16_t*)malloc(intSize*dataLength);
    if(sysEndianType==BIG_ENDIAN_SYSTEM)
    {
      memcpy(*newData, szTmpBytes+4+MetaDataByteLength+exe_params->SZ_SIZE_TYPE, dataLength*intSize);
    }
    else
    {
      unsigned char* p = szTmpBytes+4+MetaDataByteLength+exe_params->SZ_SIZE_TYPE;
      for(i=0;i<dataLength;i++,p+=intSize)
        (*newData)[i] = bytesToInt16_bigEndian(p);
    }   
  }

  else //confparams_dec->sol_ID==SZ
  {
    if (dim == 1)
      getSnapshotData_int16_1D(newData,r1,tdps, errBoundMode);
   
  } 

  free_TightDataPointStorageI2(tdps);
  if(confparams_dec->szMode!=SZ_BEST_SPEED && cmpSize!=4+sizeof(int16_t)+exe_params->SZ_SIZE_TYPE+MetaDataByteLength)
    free(szTmpBytes);
  return status;
}


void getSnapshotData_int32_1D(int32_t** data, size_t dataSeriesLength, TightDataPointStorageI* tdps, int errBoundMode)
{
  size_t i;

  if (tdps->allSameData) {
    int32_t value = bytesToInt32_bigEndian(tdps->exactDataBytes);
    *data = (int32_t*)malloc(sizeof(int32_t)*dataSeriesLength);
    for (i = 0; i < dataSeriesLength; i++)
      (*data)[i] = value;
  } else {
    decompressDataSeries_int32_1D(data, dataSeriesLength, tdps);
  }
}

int SZ_decompress_args_int32(int32_t** newData, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, unsigned char* cmpBytes, size_t cmpSize)
{
  int status = SZ_SCES;
  size_t dataLength = computeDataLength(r5,r4,r3,r2,r1);
  
  //unsigned char* tmpBytes;
  size_t targetUncompressSize = dataLength <<2; //i.e., *4
  //tmpSize must be "much" smaller than dataLength
  size_t i, tmpSize = 3+MetaDataByteLength+1+sizeof(int32_t)+exe_params->SZ_SIZE_TYPE;
  unsigned char* szTmpBytes;  
    
  if(cmpSize!=4+4+4+MetaDataByteLength && cmpSize!=4+4+8+MetaDataByteLength)
  {
    confparams_dec->losslessCompressor = 0 ; // is_lossless_compressed_data(cmpBytes, cmpSize);
    if(confparams_dec->losslessCompressor!=-1)
      confparams_dec->szMode = SZ_BEST_COMPRESSION;
    else
      confparams_dec->szMode = SZ_BEST_SPEED;   
    if(confparams_dec->szMode==SZ_BEST_SPEED)
    {
      tmpSize = cmpSize;
      szTmpBytes = cmpBytes;  
    }
    else if(confparams_dec->szMode==SZ_BEST_COMPRESSION || confparams_dec->szMode==SZ_DEFAULT_COMPRESSION)
    {
      if(targetUncompressSize<MIN_ZLIB_DEC_ALLOMEM_BYTES) //Considering the minimum size
        targetUncompressSize = MIN_ZLIB_DEC_ALLOMEM_BYTES; 
      tmpSize = sz_lossless_decompress(confparams_dec->losslessCompressor, cmpBytes, (unsigned long)cmpSize, &szTmpBytes, (unsigned long)targetUncompressSize+4+MetaDataByteLength+exe_params->SZ_SIZE_TYPE);//   (unsigned long)targetUncompressSize+8: consider the total length under lossless compression mode is actually 3+4+1+targetUncompressSize
      //szTmpBytes = (unsigned char*)malloc(sizeof(unsigned char)*tmpSize);
      //memcpy(szTmpBytes, tmpBytes, tmpSize);
      //free(tmpBytes); //release useless memory    
    }
    else
    {
   //   printf("Wrong value of confparams_dec->szMode in the double compressed bytes.\n");
      status = SZ_MERR;
      return status;
    } 
  }
  else
    szTmpBytes = cmpBytes;
  //TODO: convert szTmpBytes to data array.
  TightDataPointStorageI* tdps;
  int errBoundMode = new_TightDataPointStorageI_fromFlatBytes(&tdps, szTmpBytes, tmpSize);
  //writeByteData(tdps->typeArray, tdps->typeArray_size, "decompress-typebytes.tbt");
  int dim = computeDimension(r5,r4,r3,r2,r1); 
  int intSize = sizeof(int32_t);
  if(tdps->isLossless)
  {
    *newData = (int32_t*)malloc(intSize*dataLength);
    if(sysEndianType==BIG_ENDIAN_SYSTEM)
    {
      memcpy(*newData, szTmpBytes+4+MetaDataByteLength+exe_params->SZ_SIZE_TYPE, dataLength*intSize);
    }
    else
    {
      unsigned char* p = szTmpBytes+4+MetaDataByteLength+exe_params->SZ_SIZE_TYPE;
      for(i=0;i<dataLength;i++,p+=intSize)
        (*newData)[i] = bytesToInt32_bigEndian(p);
    }   
  }
  else if (dim == 1)
    getSnapshotData_int32_1D(newData,r1,tdps, errBoundMode);
  else
  {
  //  printf("Error: currently support only at most 4 dimensions!\n");
    status = SZ_DERR;
  }
  free_TightDataPointStorageI2(tdps);
  if(confparams_dec->szMode!=SZ_BEST_SPEED && cmpSize!=4+sizeof(int32_t)+exe_params->SZ_SIZE_TYPE+MetaDataByteLength)
    free(szTmpBytes);
  return status;
}



int SZ_decompress_args_double(double** newData, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, unsigned char* cmpBytes, 
size_t cmpSize, int compressionType, double* hist_data)
{
  int status = SZ_SCES;
  size_t dataLength = computeDataLength(r5,r4,r3,r2,r1);
  
  //unsigned char* tmpBytes;
  size_t targetUncompressSize = dataLength <<3; //i.e., *8
  //tmpSize must be "much" smaller than dataLength
  size_t i, tmpSize = 12+MetaDataByteLength_double+exe_params->SZ_SIZE_TYPE;
  unsigned char* szTmpBytes;
  if(cmpSize!=12+4+MetaDataByteLength_double && cmpSize!=12+8+MetaDataByteLength_double)
  {
    confparams_dec->losslessCompressor =0;  // is_lossless_compressed_data(cmpBytes, cmpSize);
    if(confparams_dec->szMode!=SZ_TEMPORAL_COMPRESSION)
    {
      if(confparams_dec->losslessCompressor!=-1)
        confparams_dec->szMode = SZ_BEST_COMPRESSION;
      else
        confparams_dec->szMode = SZ_BEST_SPEED;     
    }
    if(confparams_dec->szMode==SZ_BEST_SPEED)
    {
      tmpSize = cmpSize;
      szTmpBytes = cmpBytes;  
    } 
    else if(confparams_dec->szMode==SZ_BEST_COMPRESSION || confparams_dec->szMode==SZ_DEFAULT_COMPRESSION || confparams_dec->szMode==SZ_TEMPORAL_COMPRESSION)
    {
      if(targetUncompressSize<MIN_ZLIB_DEC_ALLOMEM_BYTES) //Considering the minimum size
        targetUncompressSize = MIN_ZLIB_DEC_ALLOMEM_BYTES;      
      tmpSize = sz_lossless_decompress(confparams_dec->losslessCompressor, cmpBytes, (unsigned long)cmpSize, &szTmpBytes, (unsigned long)targetUncompressSize+4+MetaDataByteLength_double+exe_params->SZ_SIZE_TYPE);      
      //szTmpBytes = (unsigned char*)malloc(sizeof(unsigned char)*tmpSize);
      //memcpy(szTmpBytes, tmpBytes, tmpSize);
      //free(tmpBytes); //release useless memory    
    }
    else
    {
    //  printf("Wrong value of confparams_dec->szMode in the double compressed bytes.\n");
      status = SZ_MERR;
      return status;
    } 
  }
  else
    szTmpBytes = cmpBytes;
    
 // confparams_dec->sol_ID = szTmpBytes[4+14]; //szTmpBytes: version(3bytes), samebyte(1byte), [14]:sol_ID=SZ or SZ_Transpose   
  //TODO: convert szTmpBytes to double array.
  TightDataPointStorageD* tdps;
  int errBoundMode = new_TightDataPointStorageD_fromFlatBytes(&tdps, szTmpBytes, tmpSize);

  int dim = computeDimension(r5,r4,r3,r2,r1);
  int doubleSize = sizeof(double);
  if(tdps->isLossless)
  {
    *newData = (double*)malloc(doubleSize*dataLength);
    if(sysEndianType==BIG_ENDIAN_SYSTEM)
    {
      memcpy(*newData, szTmpBytes+4+MetaDataByteLength_double+exe_params->SZ_SIZE_TYPE, dataLength*doubleSize);
    }
    else
    {
      unsigned char* p = szTmpBytes+4+MetaDataByteLength_double+exe_params->SZ_SIZE_TYPE;
      for(i=0;i<dataLength;i++,p+=doubleSize)
        (*newData)[i] = bytesToDouble(p);
    }   
  }
 
  else //confparams_dec->sol_ID==SZ
  {
    if(tdps->raBytes_size > 0) //v2.0
    {
      if (dim == 1)
        getSnapshotData_double_1D(newData,r1,tdps, errBoundMode, 0, hist_data);
      
      else
      {
      //  printf("Error: currently support only at most 4 dimensions!\n");
        status = SZ_DERR;
      } 
    }
    else //1.4.13 or time-based compression
    {
      if (dim == 1)
        getSnapshotData_double_1D(newData,r1,tdps, errBoundMode, compressionType, hist_data);
          
      else
      {
      //  printf("Error: currently support only at most 4 dimensions!\n");
        status = SZ_DERR;
      }     
    }
  } 

  free_TightDataPointStorageD2(tdps);
  if(confparams_dec->szMode!=SZ_BEST_SPEED && cmpSize!=12+MetaDataByteLength_double+exe_params->SZ_SIZE_TYPE)
    free(szTmpBytes); 
  return status;
}

void *SZ_decompress(int dataType, unsigned char *bytes, size_t byteLength, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
  if(confparams_dec==NULL)
    confparams_dec = (sz_params*)malloc(sizeof(sz_params));
  memset(confparams_dec, 0, sizeof(sz_params));
  if(exe_params==NULL)
    exe_params = (sz_exedata*)malloc(sizeof(sz_exedata));
  memset(exe_params, 0, sizeof(sz_exedata));
  exe_params->SZ_SIZE_TYPE = 8;

  int x = 1;
  char *y = (char*)&x;
  if(*y==1)
    sysEndianType = LITTLE_ENDIAN_SYSTEM;
  else //=0
    sysEndianType = BIG_ENDIAN_SYSTEM;

  if(dataType == SZ_FLOAT)
  {
    float *newFloatData;
    SZ_decompress_args_float(&newFloatData, r5, r4, r3, r2, r1, bytes, byteLength);
    return newFloatData;
  }
  
  else  if(dataType == SZ_INT8)
  {
    int8_t *newInt8Data;
    SZ_decompress_args_int8(&newInt8Data, r5, r4, r3, r2, r1, bytes, byteLength);
    return newInt8Data;
  }
 
  else if(dataType == SZ_INT16)
  {
    int16_t *newInt16Data;
    SZ_decompress_args_int16(&newInt16Data, r5, r4, r3, r2, r1, bytes, byteLength);
    return newInt16Data;
  }
  else if(dataType == SZ_INT32)
  {
    int32_t *newInt32Data;
    SZ_decompress_args_int32(&newInt32Data, r5, r4, r3, r2, r1, bytes, byteLength);
    return newInt32Data;
  }
  else if(dataType == SZ_DOUBLE)
  {
    double *newDoubleData;
    SZ_decompress_args_double(&newDoubleData, r5, r4, r3, r2, r1, bytes, byteLength, 0, NULL);
    return newDoubleData; 
  }
 
  else
  {
   // printf("Error: data type cannot be the types other than SZ_FLOAT or SZ_DOUBLE\n");
    return NULL;
  }
}

float rms4float(float *v, int n)
{
  int i;
  float sum = 0.0;
  for(i = 0; i < n; i++)
    sum += v[i] * v[i];
  return sqrt(sum / n);
}
float rms4double(double *v, int n)
{
  int i;
  float sum = 0.0;
  for(i = 0; i < n; i++)
    sum += v[i] * v[i];
  return sqrt(sum / n);
}
float rms4int8(int8_t *v, int n)
{
  int i;
  float sum = 0.0;
  for(i = 0; i < n; i++)
    sum += v[i] * v[i];
  return sqrt(sum / n);
}
float rms4int16(int16_t *v, int n)
{
  int i;
  float sum = 0.0;
  for(i = 0; i < n; i++)
    sum += v[i] * v[i];
  return sqrt(sum / n);
}
float rms4int32(int32_t *v, int n)
{
  int i;
  float sum = 0.0;
  for(i = 0; i < n; i++)
    sum += v[i] * v[i];
  return sqrt(sum / n);
}



#ifdef __cplusplus
}
#endif
