#include <nmod_poly_mat.h>
#include "m_basis.h"
#include "time.h"

/** Bench of polmatmul for given dimension, degree, modulus.
 *
 *  This function takes the parameters (matrix dimensions, degree bound,
 *  and prime modulus), creates two random matrices a and b, and measures
 *  the time it takes to compute c = a*b.
 *
 * \param rdim row dimension
 * \param idim inner dimension
 * \param cdim column dimension
 * \param deg degree bound
 * \param prime prime modulus
 * \param state random state for generating random matrix entries
 * \return void
 */


#define NUMBER_M_BASIS 5

static void (*m_basis[NUMBER_M_BASIS])(nmod_poly_mat_t, int64_t * ,
				       const nmod_poly_mat_t, uint64_t, const int64_t*) =
  {
    M_basis, M_basisII, M_basisIII, M_basisIV, M_basisV
  };

static char* nb_m_basis[NUMBER_M_BASIS] =  {"I", "II", "III", "IV", "V"};

void benchmark_m_basis(slong rdim, slong cdim, slong sigma, slong len,
			 ulong prime, flint_rand_t state)
{
    // create random matrices
  nmod_poly_mat_t mat;
  
  flint_randinit(state);
  nmod_poly_mat_init(mat, rdim, cdim, prime);
  nmod_poly_mat_randtest(mat, state, len + 1);
  flint_randclear(state);

  nmod_poly_mat_t res_m_basis;
  slong shifts[rdim], res_shifts[rdim];
  
  for (slong i = 0; i < rdim; i++)
    shifts[i] = 0;
    
  nmod_poly_mat_init(res_m_basis, rdim, rdim, prime);
 
  // parameters for measuring time
  double t = 0.0;
  clock_t tt;
  long nb_iter = 0;

  // let's go
  for (int i = 0; i < NUMBER_M_BASIS; i++)
    {
      t = 0.0;
      nb_iter = 0;

      while (t<0.5)
	{
	  tt = clock();
	  m_basis[i](res_m_basis, res_shifts, mat, sigma, shifts);
	  t += (double)(clock()-tt) / CLOCKS_PER_SEC;
      ++nb_iter;
	}
      t /= nb_iter;
      printf("%s\t%ld\t%ld\t%ld\t%ld\t%f\n", nb_m_basis[i], rdim, cdim, sigma, len, t);
    }
  
}

/** Launches a series of benchmarks for a given prime size.
 *
 *  Launches benchmark for many matrix dimensions and degrees for a given
 *  bitlength for the prime defining the field of coefficients.
 *
 * \param nbits bitlength of the prime modulus
 * \return void
 */
void benchmark_nbits(ulong nbits)
{
    flint_rand_t state;
    flint_randinit(state);
    const ulong prime = n_randprime(state, nbits, 0);
    flint_randclear(state);

    slong rdims[] = { 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024 };
    slong sigmas[] = { 4, 8, 16, 32 };

    printf("Bench M Basis\n");
    printf("nbits=%ld, prime=%ld\n", nbits, prime);
    printf("n\trdim\tcdim\torder\tdeg\ttime\n");
    for (size_t i = 0; i < sizeof(rdims) / sizeof(rdims[0]); ++i)
    {
        long rdim = rdims[i];
   
        for (size_t j = 0; j < sizeof(sigmas) / sizeof(sigmas[0]); ++j)
        {
            long ord = sigmas[j];
	    benchmark_m_basis(rdim, rdim / 2, ord, ord, prime, state);
	    printf("\n");
	}
    }
}

void benchmark_nbits_dim_deg(ulong nbits, ulong dim, ulong deg)
{
    flint_rand_t state;
    flint_randinit(state);
    const ulong prime = n_randprime(state, nbits, 0);
    flint_randclear(state);

    printf("Bench M Basis:\n");
    printf("nbits=%ld, dim=%ld, deg=%ld\n",nbits, dim, deg);
    printf("nbits=%ld, prime=%ld, dim=%ld, deg=%ld\n",nbits,prime,dim,deg);
    printf("n\trdim\tcdim\torder\tdeg\ttime\n");
    benchmark_m_basis(dim, dim / 2, deg, deg, prime,state);
}


int main(int argc, char *argv[])
{
    if (argc!=2 && argc!=4)
    {
        printf("Usage: %s nbits OR %s nbits rdim sigma\n",argv[0],argv[0]);
        return 1;
    }

    if (argc==2)
    {
        ulong nbits = atoi(argv[1]);
        benchmark_nbits(nbits);
        return 0;
    }

    if (argc==4)
    {
        ulong nbits = atoi(argv[1]);
        ulong dim = atoi(argv[2]);
        ulong deg = atoi(argv[3]);
        benchmark_nbits_dim_deg(nbits,dim,deg);
    }
    
    return 0;
}

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
