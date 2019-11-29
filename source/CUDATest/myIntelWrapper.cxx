#include <stdio.h>
 
void cu_vecAdd_wrapper(float*, int*);
void cu_initBuf(float*,int*);
void cu_destroy(); 
 
 extern "C" {
	void icl_vecAddTest(float*, int*);
	void icl_initBuf(float*,int*);
	void icl_destroy();
}

void icl_vecAddTest(float* val, int* n)
{	
	cu_vecAdd_wrapper( val, n );
}

void icl_initBuf( float* val, int* n )
{
	cu_initBuf(val,n);
}

void icl_destroy()
{
	cu_destroy();
}

 
 