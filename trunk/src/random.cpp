#include <stdio.h>
#include <math.h>
#include "random.h"

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

/* 
   ======================================================================  
   MT19937-2002 Algorithm

   Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


   Any feedback is very welcome.
   http://www.math.keio.ac.jp/matumoto/emt.html
   email: matumoto@math.keio.ac.jp

*/

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UMASK 0x80000000UL /* most significant w-r bits */
#define LMASK 0x7fffffffUL /* least significant r bits */
#define MIXBITS(u,v) ( ((u) & UMASK) | ((v) & LMASK) )
#define TWIST(u,v) ((MIXBITS(u,v) >> 1) ^ ((v)&1UL ? MATRIX_A : 0UL))

static uint32_t state[N]; /* the array for the state vector  */

static int leqq = 1;
static uint32_t *next;

/* CS:  From Numerical Recipes in C.  Has period of length 2.3 * 10^18 */

double ran1(long *idum) {
    int j;
    long k;
    static long iy=0;
    static long iv[NTAB];
    double temp;

    if(*idum <= 0 || !iy) {
		if(-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		for(j=NTAB+7;j>=0;j--) {
	    	k=(*idum)/IQ;
	    	*idum=IA*(*idum-k*IQ)-IR*k;
	    	if(*idum < 0) *idum+=IM;
	    	if(j<NTAB) iv[j] = *idum;
		}
		iy=iv[0];
    }
    k=(*idum)/IQ;
    *idum=IA*(*idum-k*IQ)-IR*k;
    if(*idum < 0) *idum+= IM;
    j = iy/NDIV;
    iy=iv[j];
    iv[j] = *idum;
    if((temp=AM*iy) > RNMX) return RNMX;
    else return temp;

}

void mt_srand(uint32_t s)
{
    int j;
    state[0]= s & 0xffffffffUL;
    for (j=1; j<N; j++) {
        state[j] = (1812433253UL * (state[j-1] ^ (state[j-1] >> 30)) + j); 
        //state[j] &= 0xffffffffUL;  /* for >32 bit machines */
    }
    leqq = 1;
}

// Default seed procedure
static struct mt_rand_init
{
	mt_rand_init() { mt_srand(5489UL); }
} drinit;

void mt_srand(uint32_t uKeys[], size_t uLen)
{
    uint32_t i, j;
	size_t k;
    mt_srand(19650218UL);
    i=1; j=0;
    k = (N>uLen ? N : uLen);
    for (; k; k--) {
        state[i] = (state[i] ^ ((state[i-1] ^ (state[i-1] >> 30)) * 1664525UL))
          + uKeys[j] + j; /* non linear */
        //state[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { state[0] = state[N-1]; i=1; }
        if (j>=uLen) j=0;
    }
    for (k=N-1; k; k--) {
        state[i] = (state[i] ^ ((state[i-1] ^ (state[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        //state[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { state[0] = state[N-1]; i=1; }
    }

    state[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
    leqq = 1;
}

inline void mt_next_state()
{
	leqq = N;
    next = state;
	uint32_t *p=state;
    for (int j=N-M+1; --j; p++) 
        *p = p[M] ^ TWIST(p[0], p[1]);
    for (int j=M; --j; p++) 
        *p = p[M-N] ^ TWIST(p[0], p[1]);
    *p = p[M-N] ^ TWIST(p[0], state[0]);
}

/* generates a random number on [0,0xffffffff]-interval */
uint32_t mt_rand()
{
    if (--leqq == 0)
		mt_next_state();

	uint32_t y = *next++;

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);
    return y;
}
