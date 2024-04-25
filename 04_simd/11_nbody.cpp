#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <x86intrin.h>

int main() {
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N];
  for(int i=0; i<N; i++) {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }
    __m512 x_vec, y_vec, m_vec, rx_vec, ry_vec, r_vec, r3_vec, fx_vec, fy_vec, zero_vec;
    zero_vec = _mm512_setzero_ps();
    m_vec = _mm512_load_ps(m);

    int idx[N];
    for(int i = 0;i < N;i++)idx[i] = i;
    __m512i idxvec = _mm512_loadu_si512((void*)idx);
    for(int i = 0; i < N ; i++){
	__m512 rxjvec, ryjvec, fxivec, fyivec, rsqrtvec;
	fxivec = _mm512_setzero_ps();
	fyivec = _mm512_setzero_ps();
	rxjvec = _mm512_sub_ps(
	       		_mm512_set1_ps(x[i]),
		_mm512_load_ps(x)
			);
	ryjvec = _mm512_sub_ps(
			_mm512_set1_ps(y[i]),
			_mm512_load_ps(y)
			);

	rsqrtvec = _mm512_rsqrt14_ps(
			_mm512_add_ps(
				_mm512_mul_ps(rxjvec, rxjvec),
				_mm512_mul_ps(ryjvec, ryjvec)
				)
			);

	rsqrtvec = _mm512_mul_ps(_mm512_mul_ps(rsqrtvec, rsqrtvec), rsqrtvec);

	
//	__mmask8 mask = _mm512_cmp_ps_mask(_mm512_set1_epi32(i), idxvec, _MM_CMPINT_NE);
	__mmask8 mask = _mm512_cmp_epi32_mask(idxvec, _mm512_set1_epi32(i), _MM_CMPINT_NE);
	rsqrtvec = _mm512_mask_blend_ps(mask, rsqrtvec, _mm512_setzero_ps());
	__m512 force_vec = _mm512_mul_ps(m_vec, rsqrtvec);

        //fxivec = _mm512_mask_mul_ps(zero_vec, mask, rxjvec, force_vec);
        //fyivec = _mm512_mask_mul_ps(zero_vec, mask, ryjvec, force_vec);

	fxivec =  _mm512_mul_ps(	rxjvec, force_vec);
	fyivec =  _mm512_mul_ps(        ryjvec, force_vec);
	fx[i] = -_mm512_reduce_add_ps(fxivec);
	fy[i] = -_mm512_reduce_add_ps(fyivec);
	printf("%d %g %g\n",i,fx[i],fy[i]);
    
    }

    return 0;
}

