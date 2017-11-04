#include <stdlib.h>
#include <assert.h>
#include "plcp.h"
#include <stdio.h>

void rmq_preprocess( INT * m, INT * v, INT n ) 
{
  INT i, j;
  INT lgn = flog2(n); 
  for (i = 0; i < n; i++)
    m[i*lgn] = i;
  for (j = 1; ((( INT ) 1) << j) <= n; j++)
    for (i = 0; i + (1 << j) - 1 < n; i++)
      if (v[m[i*lgn + j - 1]] < v[m[(i + (1 << (j - 1)))*lgn + j - 1]])
        m[i*lgn + j] = m[i*lgn + j - 1];
      else
        m[i*lgn + j] = m[(i + (1 << (j - 1)))*lgn + j - 1];
}
