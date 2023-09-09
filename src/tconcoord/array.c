#include <tconcoord.h>

/*============================================*/
void array_sort(int *iarr, real *darr,int size)
{
  /* sorts an integer array and a real array
     according to the value of the real array
     in increasing order
  */
  bool changed = TRUE;
  real rdum;
  int idum;
  int i;
  while(changed){
    changed = FALSE;
    for(i=0;i<size-1;i++){
      if(darr[i]>darr[i+1]){
        rdum = darr[i];
        idum = iarr[i];
        darr[i] = darr[i+1];
        iarr[i] = iarr[i+1];
        darr[i+1] = rdum;
        iarr[i+1] = idum;
        changed = TRUE;
      }
    }
  }
}
/*============================================*/
void copy_iarr(int *dest, int *src, int size)

/* copy an integer array */

{
  int i;
  for(i=0;i<size;i++){
    dest[i] = src[i];
  }
}
/*============================================*/
void i_sort(int *array, int size)
{

  bool changed = TRUE;
  int idum;
  int i;
  while(changed){
    changed = FALSE;
    for(i=0;i<size-1;i++){
      if(array[i]>array[i+1]){
        idum = array[i];
        array[i] = array[i+1];
        array[i+1] = idum;
        changed = TRUE;
      }
    }
  }
}

/*============================================================*/

void switch_real(real *x, real *y)
{
  real dummy;
  dummy=*x;
  *x=*y;
  *y=dummy;
}
/*============================================================*/

void switch_int(int *x, int *y)
{
  int dummy;
  dummy=*x;
  *x=*y;
  *y=dummy;
}

