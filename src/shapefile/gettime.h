
//Getting computation times
#ifdef WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

inline double getTime() //in millisecond
{
#ifdef WIN32
    return timeGetTime();
#else //assuming unix-type systems
    //timezone tz;
    timeval  tv;
    gettimeofday(&tv, NULL);
    return (tv.tv_sec*1000000+tv.tv_usec)*1.0/1000;
#endif
}