#include "wiedemann_algo.h"


void krylov_sequence(nmod_mat_t res, const nmod_mat_t A, const nmod_mat_t v,
		     slong n)
{
  slong rdim = A->r, cdim = A->c, vdim = v->c;
  mp_limb_t prime = A->mod.n;
  nmod_mat_t Av, A_sq, B;
  slong i, j, bound_rows;

  nmod_mat_init(Av, rdim, cdim, prime);
  nmod_mat_mul(Av, A, v);

  for (i = 0; i < rdim; i++)
    {
      _nmod_vec_set(res->rows[i], A->rows[i], vdim);
      _nmod_vec_set(res->rows[i + rdim], Av->rows[i], vdim);
    }

  nmod_mat_clear(Av);
  
  for (i = 0; i < n; i++)
    {
      bound_rows = ((slong) pow(2, i + 1)) * rdim;
      nmod_mat_mul(A_sq, A, A);
      nmod_mat_window_init(Av, res, 0, 0, bound_rows - 1 , vdim);
      nmod_mat_init_set(B, Av);
      nmod_mat_mul(B, A_sq);
      for (j = 0; j < B->r; j++)
	_nmod_vec_set(res->rows[bound_rows + j], Av->rows[j], vdim);

      nmod_mat_window_clear(Av);
      nmod_mat_clear(B);
    }
}

void Berlekamp_Massey(nmod_poly_t res, const nmod_poly_t P);

void Wiedemann(nmod_poly_t res, const nmod_mat_t A);

void Block_Wiedemann(nmod_poly_t res, const nmod_mat_t A);
