/* Linux hostnm routine */
 
#include <unistd.h>
 
/* function names link differently depending on OS */
 
#include "fpreproc/f77name.h"
 
#ifdef __LINUX
int F77NAME(hostnm)(name,len)
#else
int F77NAME(hostnm_dummy)(name,len)
#endif
char *name;
int len;
{
      return gethostname(name,len);
}
