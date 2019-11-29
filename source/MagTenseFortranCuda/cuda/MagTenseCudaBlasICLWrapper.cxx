#include <stdio.h>
 
void cu_initDemagMatrices( const float*, const float*, const float*, const float*, const float*, const float*, int*);
void cu_MVMult_GetH( const float*, const float*, const float*, float*, float*, float*, int*, float* );
void cu_destroy();
void cu_test( const float*, int* n);

 extern "C" {
	void icl_initDemagMatrices(const float*, const float*, const float*, const float*, const float*, const float*, int*);
	
	void icl_MVMult_GetH( const float*, const float*, const float*, float*, float*, float*, int*, float* );
	
	void icl_destroy();
	
	void icl_test( const float* K, int* n);
}

void icl_initDemagMatrices(const float* Kxx, const float* Kxy, const float* Kxz, const float* Kyy, const float* Kyz, const float* Kzz, int* n)
{	
	cu_initDemagMatrices( Kxx, Kxy, Kxz, Kyy, Kyz, Kzz, n );
}

void icl_MVMult_GetH( const float* Mx, const float* My, const float* Mz, float* Hx, float* Hy, float* Hz, int* n, float* pref )
{
	cu_MVMult_GetH( Mx, My, Mz, Hx, Hy, Hz, n, pref );
}

void icl_destroy()
{
	cu_destroy();
}

 void icl_test(const float* K, int* n)
 {
	 cu_test(K, n);
 }
 