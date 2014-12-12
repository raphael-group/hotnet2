#include <math.h>

#define index_in(a, i, j) a[i*p + j]
#define index_out(a, i, j) a[i*q + j]

void compute_sim(double *infmat, double *h, long *indices, long p, long q, double *M)
{
  int i, j;
  for (i = 0; i < q; ++i)
  {
    for (j = 0; j < q; ++j)
    {
      index_out(M, i, j) = index_in(infmat, indices[i], indices[j]) * h[j];
    }
  }
}

void compute_sim_classic(double *infmat, double *h, long *indices, long p, long q, double *M)
{
  int i, j;
  for (i = 0; i < q; ++i)
  {
    for (j = 0; j < q; ++j)
    {
      index_out(M, i, j) = index_out(M, j, i) = fmin(index_in(infmat, indices[i], indices[j]),
                                                     index_in(infmat, indices[j], indices[i])) 
                                                * fmax(h[i], h[j]);
    }
  }
}
