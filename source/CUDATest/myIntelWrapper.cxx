#include <stdio.h>
 
void vecAdd_wrapper(void);
 
 extern "C" {
	void vecAddTest(void);
}

void vecAddTest(void)
{
	vecAdd_wrapper();
}
 
 