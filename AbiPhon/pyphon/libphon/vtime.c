/*****************************************************************
      FORTRAN interface to get user-time and real time on
      AIX-system might work on other UNIX systems as well
      (must work on all BSDish systems and SYSVish systems
      with BSD-compatibility timing routines ...)!?
*****************************************************************/

#ifndef NULL
#define NULL    ((void *)0)
#endif

#include <sys/times.h>
#include <sys/resource.h>

#include <portinfo>

int    gettimeofday();
int    getrusage();

#if defined (F77EXTERNS_COMPAQ_F90)
#define vtime vtime_
#elif defined (F77EXTERNS_EXTRATRAILINGBAR)
#define vtime vtime__
#elif defined (F77EXTERNS_LOWERCASE_TRAILINGBAR)
#define vtime vtime_
#elif defined (F77EXTERNS_NOTRAILINGBAR)
#define vtime vtime
#elif defined (F77EXTERNS_UPPERCASE_NOTRAILINGBAR)
#define vtime VTIME
#endif

void vtime(double *cputim, double *wall)
{
        struct rusage ppt;
        struct rusage cpt;
        struct timeval tpu,tps;
        struct timeval tcu,tcs;
        struct timeval now;
        int    ierr;
        ierr = getrusage(RUSAGE_SELF,&ppt);
        tpu  = ppt.ru_utime;
        tps  = ppt.ru_stime;
        ierr = getrusage(RUSAGE_CHILDREN,&cpt);
        tcu  = cpt.ru_utime;
        tcs  = cpt.ru_stime;
        ierr = gettimeofday(&now,NULL);
        *cputim=((double) tpu.tv_sec) + ((double) tpu.tv_usec) / 1e6 +
                ((double) tps.tv_sec) + ((double) tps.tv_usec) / 1e6 +
                ((double) tcu.tv_sec) + ((double) tcu.tv_usec) / 1e6 +
                ((double) tcs.tv_sec) + ((double) tcs.tv_usec) / 1e6;
        *wall  =((double) now.tv_sec) + ((double) now.tv_usec) / 1e6;
}

