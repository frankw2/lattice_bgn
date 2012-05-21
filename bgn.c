#include "bgn.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cmath>

// prime represents the prime we are going to mod by for ciphertext
ZZ prime;

//prime is 2^k
long k;

//number of total records in the database (how big the database is)
ZZ num_Records; 

//the n,m parameter in the paper
long n,m;

//the public and private keys
mat_ZZ a,t;

void GenerateKeys()
{
	//generate a and t
  // Primitive matrix G
  mat_ZZ G, A_bar, H, R, A, temp;
  long w = n*k;
  G.SetDims(n,w);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < k; j++) {
      G[i][j] = (int) pow(2.0f, j);
    }
  }

  long m_bar = 2*n;
  A_bar.SetDims(n,m_bar);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i == j) {
        A_bar[i][j] = 1;
      } else {
        A_bar[i][j] = 0;
      }
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = n+1; j < m_bar; j++) {
      A_bar[i][j] = RandomBnd(prime);
    }
  }
  
  H = ident_mat_ZZ(n);

  const gsl_rng_type * T;
  gsl_rng * r;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  R.SetDims(m_bar,w);
  for (int i = 0; i < m_bar; i++) {
    for (int j = 0; j < w; j++) {
      R[i][j] = (int) gsl_ran_gaussian(r, n);
      R[i][j] %= prime;
    }
  }

  m = m_bar + w;
  A.SetDims(n,m);
	temp = H*G - A_bar*R;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m_bar; j++) {
      A[i][j] = A_bar[i][j];
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < w; j++) {
      A[i][m_bar + j] = temp[i][j];
    }
  }
	//also set parameters based on restrictions in the paper (should just do set values more or less?)
}

//message should be square matrix m x m 
void Encrypt(mat_ZZ& ciphertext, mat_ZZ& message)
{
	//set and create random matrix
	mat_ZZ s;
	s.SetDims(n,m);
	for (int i = 0; i < n; i++){
		for (int j = 0; j < m; j++){
			s[i][j] = RandomBnd(prime);
		}
	}
	
	mat_ZZ term1;
	term1 = a*s;
	
	ciphertext = term1 + message;
	
	//need to zero out least significant bits (don't include num_Records significant bits)
	
}

void Decrypt(mat_ZZ& message, mat_ZZ& ciphertext)
{
	mat_ZZ e;
	
	//set dimensions for e?
	
	e = t*ciphertext*transpose(t);
	//mod all elements by the prime
	for (int i = 0; i < e.NumRows(); i++){
		for (int j = 0; j < e.NumCols(); j++){
			e[i][j] = e[i][j] % prime;
		}
	}
	message = inv(t)*e*inv(transpose(t));
	for (int i = 0; i < message.NumRows(); i++){
		for (int j = 0; j < message.NumCols(); j++){
			message[i][j] = message[i][j] % num_Records;
		}
	}
}
