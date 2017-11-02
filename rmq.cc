#include <stdlib.h>
#include <assert.h>
#include "plcp.h"
#include <stdio.h>


// $m$ must have size $n\lg n$ and $a$ is a vecter of size $n$
void rmq_preprocess(INT *m, INT *v, INT n) {

  INT i, j;
  INT lgn = flog2(n); 
   
  // initialize $m$ for the intervals with length $1$
  for (i = 0; i < n; i++)
    m[i*lgn] = i;
  // compute values from smaller to bigger intervals
  for (j = 1; ((( INT ) 1) << j) <= n; j++)
    for (i = 0; i + (1 << j) - 1 < n; i++)
      if (v[m[i*lgn + j - 1]] < v[m[(i + (1 << (j - 1)))*lgn + j - 1]])
        m[i*lgn + j] = m[i*lgn + j - 1];
      else
        m[i*lgn + j] = m[(i + (1 << (j - 1)))*lgn + j - 1];
}

INT rmq(INT *m, INT *v, INT n, INT i, INT j) {
  INT lgn = flog2(n); 
  if (i > j) {INT tmp = j; j = i; i = tmp;}
  i++; //for LCP
  if (i == j) return i;
  INT k = flog2(j-i+1); 
  INT a = m[i * lgn + k]; 
  INT shift = ((( INT ) 1) << k);
  INT b = m[(j - shift + 1) * lgn + k];  
  return v[a]>v[b]?b:a;
}

#if 0

#include <stdlib.h>

#include <assert.h>

INT main() {
  INT a[] = {2, 1, 1, 4, 0, 5, 3, -1, 3, 3, 2};
  INT n = 10, lgn = flog2(n);

  INT *m = malloc(sizeof(INT)*n*lgn);
  rmq_preprocess(m, a, n);

  INT x = 0;
  INT y = 3;

  printf("%d\n", rmq(m, a, n, x, y));

  free(m);

  return EXIT_SUCCESS;
}

#endif
