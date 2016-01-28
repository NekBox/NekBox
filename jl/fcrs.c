#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "c99.h"
#include "name.h"
#include "fail.h"
#include "types.h"
#include "mem.h"
#include "comm.h"
#include "crs.h"

static struct crs_data **handle_array = 0;
static int handle_max = 0;
static int handle_n = 0;


/*--------------------------------------------------------------------------
   FORTRAN Interface to coarse solver
  --------------------------------------------------------------------------*/

#undef crs_setup_amg
#undef crs_solve_amg
#undef crs_stats_amg
#undef crs_free_amg
#define ccrs_setup   PREFIXED_NAME(crs_setup_amg)
#define ccrs_solve   PREFIXED_NAME(crs_solve_amg)
#define ccrs_stats   PREFIXED_NAME(crs_stats_amg)
#define ccrs_free    PREFIXED_NAME(crs_free_amg )

#define fcrs_setup   FORTRAN_NAME(crs_setup_amg,CRS_SETUP_AMG)
#define fcrs_solve   FORTRAN_NAME(crs_solve_amg,CRS_SOLVE_AMG)
#define fcrs_stats   FORTRAN_NAME(crs_stats_amg,CRS_STATS_AMG)
#define fcrs_free    FORTRAN_NAME(crs_free_amg ,CRS_FREE_AMG)

void fcrs_setup(sint *handle, const MPI_Fint *comm, const sint *np,
                const sint *n, const slong id[], const sint *nz,
                const sint Ai[], const sint Aj[], const double A[],
                const sint *null_space)
{
  struct comm c;
  if(handle_n==handle_max)
    handle_max+=handle_max/2+1,
    handle_array=trealloc(struct crs_data*,handle_array,handle_max);
  comm_init_check(&c, *comm, *np);
  handle_array[handle_n]=ccrs_setup(*n,(const ulong*)id,
                                    *nz,(const uint*)Ai,(const uint*)Aj,A,
                                    *null_space,&c);
  comm_free(&c);
  *handle = handle_n++;
}

#define CHECK_HANDLE(func) do \
  if(*handle<0 || *handle>=handle_n || !handle_array[*handle]) \
    fail(1,__FILE__,__LINE__,func ": invalid handle"); \
while(0)

void fcrs_solve(const sint *handle, double x[], double b[])
{
  CHECK_HANDLE("crs_solve");
  ccrs_solve(x,handle_array[*handle],b);
}

void fcrs_stats(const sint *handle)
{
  CHECK_HANDLE("crs_stats");
  ccrs_stats(handle_array[*handle]);
}

void fcrs_free(sint *handle)
{
  CHECK_HANDLE("crs_free");
  ccrs_free(handle_array[*handle]);
  handle_array[*handle] = 0;
}

// XXT version

#undef crs_setup_xxt
#undef crs_solve_xxt
#undef crs_stats_xxt
#undef crs_free_xxt
#define ccrs_setup   PREFIXED_NAME(crs_setup_xxt)
#define ccrs_solve   PREFIXED_NAME(crs_solve_xxt)
#define ccrs_stats   PREFIXED_NAME(crs_stats_xxt)
#define ccrs_free    PREFIXED_NAME(crs_free_xxt )

#define fcrs_setup   FORTRAN_NAME(crs_setup_xxt,CRS_SETUP_XXT)
#define fcrs_solve   FORTRAN_NAME(crs_solve_xxt,CRS_SOLVE_XXT)
#define fcrs_stats   FORTRAN_NAME(crs_stats_xxt,CRS_STATS_XXT)
#define fcrs_free    FORTRAN_NAME(crs_free_xxt ,CRS_FREE_XXT)


void fcrs_setup(sint *handle, const MPI_Fint *comm, const sint *np,
                const sint *n, const slong id[], const sint *nz,
                const sint Ai[], const sint Aj[], const double A[],
                const sint *null_space)
{
  struct comm c;
  if(handle_n==handle_max)
    handle_max+=handle_max/2+1,
    handle_array=trealloc(struct crs_data*,handle_array,handle_max);
  comm_init_check(&c, *comm, *np);
  handle_array[handle_n]=ccrs_setup(*n,(const ulong*)id,
                                    *nz,(const uint*)Ai,(const uint*)Aj,A,
                                    *null_space,&c);
  comm_free(&c);
  *handle = handle_n++;
}

#define CHECK_HANDLE(func) do \
  if(*handle<0 || *handle>=handle_n || !handle_array[*handle]) \
    fail(1,__FILE__,__LINE__,func ": invalid handle"); \
while(0)

void fcrs_solve(const sint *handle, double x[], double b[])
{
  CHECK_HANDLE("crs_solve");
  ccrs_solve(x,handle_array[*handle],b);
}

void fcrs_stats(const sint *handle)
{
  CHECK_HANDLE("crs_stats");
  ccrs_stats(handle_array[*handle]);
}

void fcrs_free(sint *handle)
{
  CHECK_HANDLE("crs_free");
  ccrs_free(handle_array[*handle]);
  handle_array[*handle] = 0;
}


