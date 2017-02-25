/* Linux chdir routine */
 
#include "fpreproc/f77name.h"
 
/* function names link differently depending on OS */
 
int F77NAME(linux_chdir)(name)
char *name;
{
      return chdir(name);
}
 
 
 
 
