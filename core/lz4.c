#include <stdio.h>
#include <string.h>
#include "./lz4/lib/lz4.h"

static LZ4_stream_t lz4Stream_body; 
static LZ4_stream_t* lz4Stream=&lz4Stream_body;
static LZ4_streamDecode_t lz4StreamDecode_body;
static LZ4_streamDecode_t* lz4StreamDecode = &lz4StreamDecode_body;
enum {
  BLOCK_BYTES = 8*1024
};
static char inpBuf[2][BLOCK_BYTES];
static int inpBufIndex=0;
static char decBuf[2][BLOCK_BYTES];
static int  decBufIndex = 0;

#ifndef FNAME_H
#define FNAME_H

/*
   FORTRAN naming convention
     default      cpgs_setup, etc.
     -DUPCASE     CPGS_SETUP, etc.
     -DUNDERSCORE cpgs_setup_, etc.
*/

#ifdef UPCASE
#  define FORTRAN_NAME(low,up) up
#else
#ifdef UNDERSCORE
#  define FORTRAN_NAME(low,up) low##_
#else
#  define FORTRAN_NAME(low,up) low
#endif
#endif

#endif

#define lz4_pack  FORTRAN_NAME(lz4_pack,  LZ4_PACK)
#define lz4_unpack  FORTRAN_NAME(lz4_unpack,  LZ4_UNPACK)

/** 
 * Compress input buffer in to output buffer out. 
 * The input and output buffer are assumed allocated and of size count with each
 * element being of size n.
 *
 * @param in input buffer
 * @param out output buffer
 * @param n number of buffer elements
 * @param size size of each elmement
 * @param count number of elements
 * @param ierr error code
 */

void lz4_pack(const void * const in, const int * const size, void * const out, int * const compSize, const int * const ierr) {
  int i=0;

  size_t ccount=0; // Bytes in the chunk to be compressed. This is BLOCK_BYTES besides for the last chunk
  LZ4_resetStream(lz4Stream);
  if(*size<0 || out==NULL || *size<0 ) {
    printf("Error in Nek LZ4 compression\n");
    return;
  }
  /*printf("size: %d\n", *size);*/
  size_t lcount=(size_t) *size; // Bytes left to be compressed
  /*printf("lcount: %zu\n", lcount);*/
  size_t offset=0;
  *compSize=0;
  while (lcount > 0) {
    if(lcount < BLOCK_BYTES) {
      ccount=lcount;
    }
    else {
      ccount=BLOCK_BYTES;
    }
    /*printf("offset: %zu\n",offset);*/
    /*printf("ccount: %zu\n",ccount);*/
    char* const inpPtr = inpBuf[inpBufIndex];
    memcpy(inpPtr, in+offset, ccount);
    /*printf("i: %d\n",i);*/
    {
      char cmpBuf[LZ4_COMPRESSBOUND(BLOCK_BYTES)];
      const int cmpBytes = LZ4_compress_fast_continue(
          lz4Stream, inpPtr, cmpBuf, ccount, sizeof(cmpBuf), 1);
      if(cmpBytes <= 0) {
        break;
      }
      /*printf("cmpBytes: %d\n",cmpBytes);*/
      /*printf("compSize: %d\n",*compSize);*/
      memcpy(out+*compSize, &cmpBytes, sizeof(cmpBytes));
      memcpy(out+*compSize+sizeof(cmpBytes), &cmpBuf, cmpBytes);
      *compSize=*compSize+cmpBytes+sizeof(cmpBytes);
    }
    inpBufIndex = (inpBufIndex + 1) % 2;
    lcount=lcount-ccount;
    offset=offset+ccount;
    /*printf("----------------------\n");*/
    i=i+1;
  }
  /*printf("compSize in C: %d\n", *compSize);*/
}

void lz4_unpack(void * in, const size_t * const compSize, void * const out, const int * const size, const int * const ierr) {
  size_t offset=0;
  size_t offset_in=0;
  decBufIndex=0;
  LZ4_setStreamDecode(lz4StreamDecode,NULL,0);
  for(;;) {
    char cmpBuf[LZ4_COMPRESSBOUND(BLOCK_BYTES)];
    int cmpBytes = 0;
    /*printf("offset_in+sizeof(cmpBytes): %d\n",offset_in+sizeof(cmpBytes));*/
    /*printf("compSize: %d\n",*compSize);*/
    int tmp=offset_in+sizeof(cmpBytes);
    /*printf("tmp: %d\n", tmp);*/
    /*printf("compSize: %d\n", *compSize);*/
    if(tmp > (int) *compSize) {
      break;
    }
    memcpy(&cmpBytes, in+offset_in, sizeof(cmpBytes));
    /*printf("cmpBytes: %d\n",cmpBytes);*/
    if(cmpBytes <= 0) {
      /*printf("cmpBytes %d\n", cmpBytes);*/
      break;
    }

    /*printf("lz4 cmpBytes %d\n", cmpBytes);*/
    memcpy(&cmpBuf, in+offset_in+sizeof(cmpBytes), cmpBytes);
    const size_t readCount1=0;
    offset_in=offset_in+sizeof(cmpBytes)+cmpBytes;
    /*printf("new offset_in: %zu\n", offset_in);*/

    {
      char* decPtr = decBuf[decBufIndex];
      int decBytes = LZ4_decompress_safe_continue(
          lz4StreamDecode, cmpBuf, decPtr, cmpBytes, BLOCK_BYTES);
      if(decBytes <= 0) {
        /*printf("lz4 decBytes: %d\n", decBytes);*/
        break;
      }
      memcpy(out+offset, decPtr, (size_t) decBytes);
      offset=offset+(size_t) decBytes;
    }
    decBufIndex = (decBufIndex + 1) % 2;
  }
}
