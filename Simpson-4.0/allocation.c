/*
 * allocation.c
 *
 *  Created on: 22.5.2012
 *      Author: zdenek
 */
#include <stdio.h>
#include <stdlib.h>
#include <tcl.h>

#define MEMALLOCATOR(VAR,TYPE,SIZE,ERRORMESSAGE) \
  TYPE *VAR = ( TYPE *)Tcl_Alloc( (SIZE) * sizeof(TYPE) ); \
	if ( VAR == NULL ) { \
		fprintf(stderr,"Error in allocation: %s\n",ERRORMESSAGE); \
		exit(1); \
	}
#define MEMRELEASE(VAR) \
		Tcl_Free((char*)VAR);


void test(void)
{

MEMALLOCATOR(a,int,5,"kkk");

MEMRELEASE(a);


}
