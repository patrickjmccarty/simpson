/*
 * timing.h
 *
 *  Created on: 14.6.2012
 *      Author: Zdenek
 *
 *  Timing macros for windows
 */

#ifdef TIMING

#include <windows.h>
#include <winbase.h>


#define TIMING_INIT \
	LARGE_INTEGER _tickpsec_;\
	QueryPerformanceFrequency(&_tickpsec_);

#define TIMING_INIT_VAR(VAR)\
	LARGE_INTEGER VAR;
#define TIMING_INIT_VAR2(v1,v2)\
		LARGE_INTEGER v1, v2, _tickpsec_;\
		QueryPerformanceFrequency(&_tickpsec_);

#define TIMING_TIC(v)\
	QueryPerformanceCounter(&v);

#define TIMING_TOC(tv1,tv2,MESSAGE)\
		QueryPerformanceCounter(&tv2);\
		printf("Timing %s: %.9f\n",MESSAGE,((float)(tv2.QuadPart-tv1.QuadPart))/(float)_tickpsec_.QuadPart);

#else

#define TIMING_INIT
#define TIMING_INIT_VAR(v)
#define TIMING_INIT_VAR2(v1,v2)
#define TIMING_TIC(v)
#define TIMING_TOC(v1,v2,m)

#endif
