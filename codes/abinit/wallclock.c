#include <time.h>
#include <sys/time.h>

//To use in fortran
//  INTERFACE 
//     FUNCTION wallclock()
//       integer, parameter :: dp=kind(1.0d0)
//       REAL(dp) :: wallclock
//     END FUNCTION wallclock
//  END INTERFACE
// ...
// ...
// integer, parameter :: dp=kind(1.0d0)
// real(dp) :: t1,t2
// t1 = wallclock()
//    call myfct(...)
// t2 = wallclock()

//For FORTRAN
double wallclock_(){
  struct timeval timer;
  gettimeofday(&timer, NULL);
  double time = timer.tv_sec + timer.tv_usec * 1.0E-6;
  return time;
}


//For C
double wallclock(){
  struct timeval timer;
  gettimeofday(&timer, NULL);
  double time = timer.tv_sec + timer.tv_usec * 1.0E-6;
  return time;
}

