#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include <sys/types.h>
#ifdef LZ4COMPRESSION
#include "./lz4/lib/lz4.h"
#endif

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

#define byte_reverse  FORTRAN_NAME(byte_reverse,  BYTE_REVERSE)
#define byte_reverse8 FORTRAN_NAME(byte_reverse8, BYTE_REVERSE8)
#define byte_open     FORTRAN_NAME(byte_open,     BYTE_OPEN   )
#define byte_close    FORTRAN_NAME(byte_close,    BYTE_CLOSE  )
#define byte_rewind   FORTRAN_NAME(byte_rewind,   BYTE_REWIND )
#define byte_read     FORTRAN_NAME(byte_read,     BYTE_READ   )
#define byte_write    FORTRAN_NAME(byte_write,    BYTE_WRITE  )
#define set_bytesw_write FORTRAN_NAME(set_bytesw_write,SET_BYTESW_WRITE)
#define set_bytesw_read  FORTRAN_NAME(set_bytesw_read ,SET_BYTESW_READ )
#define get_bytesw_write FORTRAN_NAME(get_bytesw_write,GET_BYTESW_WRITE)
#define get_bytesw_read  FORTRAN_NAME(get_bytesw_read ,GET_BYTESW_READ )

#define READ     1
#define WRITE    2
#define MAX_NAME 132

#define SWAP(a,b)       temp=(a); (a)=(b); (b)=temp;

static FILE *fp=NULL;
static int  flag=0;
static char name[MAX_NAME+1];

int bytesw_write=0;
int bytesw_read=0;

#ifdef LZ4COMPRESSION
static LZ4_stream_t lz4Stream_body; 
static LZ4_stream_t* lz4Stream=&lz4Stream_body;
static LZ4_streamDecode_t lz4StreamDecode_body;
static LZ4_streamDecode_t* lz4StreamDecode = &lz4StreamDecode_body;
enum {
  BLOCK_BYTES = 1024 * 8
};
static char inpBuf[2][BLOCK_BYTES];
static int inpBufIndex=0;
static char decBuf[2][BLOCK_BYTES];
static int  decBufIndex = 0;
static void *decStream = NULL; 
#endif
/*************************************byte.c***********************************/

#ifdef UNDERSCORE
  void exitt_();
#else
  void exitt();
#endif

#ifdef LZ4COMPRESSION
size_t fwrite_lz4 (const void *buf, size_t size, size_t count, FILE *fp) {
  int i=0;
  size_t ccount=0; // Bytes in the chunk to be compressed. This is BLOCK_BYTES besides for the last chunk
  if(count<0 || fp==NULL) {
    return -1;
  }
  size_t lcount=count*size; // Bytes left to be compressed
  size_t offset=0;
  while (lcount > 0) {
    if(lcount < BLOCK_BYTES) {
      ccount=lcount;
    }
    else {
      ccount=BLOCK_BYTES;
    }
    char* const inpPtr = inpBuf[inpBufIndex];
    memcpy(inpPtr, buf+offset, ccount);
    {
      char cmpBuf[LZ4_COMPRESSBOUND(BLOCK_BYTES)];
      const int cmpBytes = LZ4_compress_fast_continue(
          lz4Stream, inpPtr, cmpBuf, ccount, sizeof(cmpBuf), 1);
      if(cmpBytes <= 0) {
        break;
      }
      fwrite(&cmpBytes, sizeof(cmpBytes),1,fp);
      fwrite(&cmpBuf, 1, cmpBytes, fp);
    }
    inpBufIndex = (inpBufIndex + 1) % 2;
    lcount=lcount-ccount;
    offset=offset+ccount;
  }
  return count;
}

FILE* fopen_lz4 (FILE *fp) {
  size_t offset=0;
  int k=1;
  decStream=malloc(BLOCK_BYTES*100);
  for(;;) {
    char cmpBuf[LZ4_COMPRESSBOUND(BLOCK_BYTES)];
    int cmpBytes = 0;
    const size_t readCount0 = fread(&cmpBytes, sizeof(cmpBytes), 1, fp);
    if(readCount0 != 1 || cmpBytes <= 0) {
      printf("lz4 readCount0: %zu, cmpBytes %d\n", readCount0, cmpBytes);
      break;
    }

      printf("lz4 cmpBytes %d\n", cmpBytes);
    const size_t readCount1 = fread(cmpBuf,1,cmpBytes,fp);
    if(readCount1 != cmpBytes) {
      printf("lz4 readCount1: %zu, cmpBytes %d\n", readCount1, cmpBytes);
      break;
    }

    {
      char* decPtr = decBuf[decBufIndex];
      int decBytes = LZ4_decompress_safe_continue(
          lz4StreamDecode, cmpBuf, decPtr, cmpBytes, BLOCK_BYTES);
      if(decBytes <= 0) {
        printf("lz4 decBytes: %d\n", decBytes);
        break;
      }
      size_t required=offset+(size_t) decBytes;
      while(required>BLOCK_BYTES*100*k) {
        k=k+1;
        decStream=realloc(decStream,BLOCK_BYTES*100*k);
      }
      memcpy(decStream+offset, decPtr, (size_t) decBytes);
      offset=offset+(size_t) decBytes;
    }
    decBufIndex = (decBufIndex + 1) % 2;
    printf("lz4 k: %d\n", k);
  }
  fclose(fp);
  return fmemopen(decStream, BLOCK_BYTES*100*k,"rw");
}
#endif

void byte_reverse(float *buf, int *nn,int *ierr)
{
  int n;
  char temp, *ptr;

  if (*nn<0)
  {
    printf("byte_reverse() :: n must be positive\n"); 
    *ierr=1;
    return;
  }
  
  for (ptr=(char *)buf,n=*nn; n--; ptr+=4)
  {
     SWAP(ptr[0],ptr[3])
     SWAP(ptr[1],ptr[2])
  }
  *ierr=0;
}

void byte_reverse8(float *buf, int *nn,int *ierr)
{
  int n;
  char temp, *ptr;

  if (*nn<0)
  {
    printf("byte_reverse8() :: n must be positive\n");
    *ierr=1;
    return;
  }
  if(*nn % 2 != 0)
  {
    printf("byte_reverse8() :: n must be multiple of 2\n");
    *ierr=1;
    return;
  }

  for (ptr=(char *)buf,n=*nn,n=n+2; n-=2; ptr+=8)
  {
     SWAP(ptr[0],ptr[7])
     SWAP(ptr[1],ptr[6])
     SWAP(ptr[2],ptr[5])
     SWAP(ptr[3],ptr[4])
  }
  *ierr=0;
}


void byte_open(char *n,int *ierr)
{
  int  i,len,istat;
  char slash;
  char dirname[MAX_NAME+1];

  len = strlen(n);
  
  if (len<0)
  {
    printf("byte_open() :: file name has negative length!\n"); 
    *ierr=1;
    return;
  }

  if (len>MAX_NAME)
  {
    printf("byte_open() :: file name too long!\n"); 
    *ierr=1;
    return;
  }

  strcpy(name,n);
  strcpy(dirname,n);

  for (i=1;dirname[i]!='\0';i++)
  {
     if (i>0 && dirname[i]=='/')
     {
       slash = name[i];
       dirname[i] = '\0';
       istat = mkdir(dirname,0755);
     }
  }
  *ierr=0;
}

void byte_close(int *ierr)
{
  if (!fp) return;

  if (fclose(fp))
  {
    printf("byte_close() :: couldn't fclose file!\n");
    *ierr=1;
    return;
  }
#ifdef LZ4COMPRESSION
  else if(flag==READ) {
    free(decStream);
  }
#endif

  fp=NULL;
  *ierr=0;
}

void byte_rewind()
{
  if (!fp) return;

  rewind(fp);
}


void byte_write(float *buf, int *n,int *ierr)
{
  int flags;
  mode_t mode;

  if (*n<0)
  {
    printf("byte_write() :: n must be positive\n"); 
    *ierr=1;
    return;
  }

  if (!fp)
  {
    if (!(fp=fopen(name,"wb")))
    {
      printf("byte_write() :: fopen failure!\n"); 
      *ierr=1;
      return;
    }
    flag=WRITE;
#ifdef LZ4COMPRESSION
    inpBufIndex=0;
    LZ4_resetStream(lz4Stream);
#endif
  }

  if (flag==WRITE)
    {
      if (bytesw_write == 1)
        byte_reverse (buf,n,ierr);
#ifdef LZ4COMPRESSION
      fwrite_lz4(buf,sizeof(float),*n,fp);
#else
      fwrite(buf,sizeof(float),*n,fp);
#endif
    }
  else
  {
      printf("byte_write() :: can't fwrite after freading!\n"); 
      *ierr=1;
      return;
  }
  *ierr=0;
}


void byte_read(float *buf, int *n,int *ierr)
{
  int flags;
  mode_t mode;
#ifdef LZ4COMPRESSION
  char key[] = ".re2";
#endif

  if (*n<0)
    {printf("byte_read() :: n must be positive\n"); *ierr=1; return;}

  if (!fp)
  {
     if (!(fp=fopen(name,"rb")))
     {
        printf("%s\n",name);
        printf("byte_read() :: fopen failure2!\n"); 
        *ierr=1;
        return;
     }
     flag=READ;
#ifdef LZ4COMPRESSION
     decBufIndex=0;
     LZ4_setStreamDecode(lz4StreamDecode,NULL,0);
     if(strstr(name,key) == 0) {
       fp=fopen_lz4(fp);
     }
#endif
  }

  if (flag==READ)
  {
     if (bytesw_read == 1)
        byte_reverse (buf,n,ierr);
     fread(buf,sizeof(float),*n,fp);
     if (ferror(fp))
     {
       printf("ABORT: Error reading %s\n",name);
       *ierr=1;
       return;
     }
     else if (feof(fp))
     {
       printf("ABORT: EOF found while reading %s\n",name);
       *ierr=1;
       return;
     }

  }
  else
  {
     printf("byte_read() :: can't fread after fwriting!\n"); 
     *ierr=1;
     return;
  }
  *ierr=0;
}

void set_bytesw_write (int *pa)
{
    if (*pa != 0)
       bytesw_write = 1;
    else
       bytesw_write = 0;
}

void set_bytesw_read (int *pa)
{
    if (*pa != 0)
       bytesw_read = 1;
    else
       bytesw_read = 0;
}

void get_bytesw_write (int *pa)
{
    *pa = bytesw_write;
}

void get_bytesw_read (int *pa)
{
    *pa = bytesw_read;
}
