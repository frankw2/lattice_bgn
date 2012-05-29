#include "bgn.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <NTL/LLL.h>
#include <cmath>
#include <ctime>

// prime represents the prime we are going to mod by for ciphertext
ZZ prime, q;

//q is 2^k
int k;

//number of total records in the database (how big the database is)
ZZ num_Records; 

//the n,m parameter in the paper
long n,m;

//the public and private keys
mat_ZZ a,t;

const gsl_rng_type * T;
gsl_rng *r;

void Init() {
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_set(r, time(NULL));
  k = 5;
  q = power2_ZZ(k);
  n = 5;
  printf("done with init\n");
}

// Samples z such that Gz = v. G is n x nk.
void SamplePreImageG(vec_ZZ & z, mat_ZZ & G, vec_ZZ & g, int n, int k, vec_ZZ v) {
  ZZ u;
  for (int i = 0; i < n; i++) {
    u = v[i];
    vec_ZZ zpart;
    zpart.SetLength(k);
    for (int j = 0; j < k; j++) {
      int temp = (2*(int)(gsl_ran_gaussian(r, pow(log(n),0.5))) + u%2);
      u = (u-temp)/2;
      z[i*k+j] = temp;
      zpart[j] = temp;
    }
    /*cout << "======== u =======" << endl;
    cout << v[i] % q << endl;
    cout << "======== g*z =======" << endl;
    cout << g*zpart % q << endl;*/
  }
}

void SamplePreImageA(vec_ZZ & x, mat_ZZ & G, mat_ZZ & A_bar, mat_ZZ & R, mat_ZZ & RI, vec_ZZ & g, int n, int k, int m_bar, int w, vec_ZZ u) {
  vec_ZZ p, p1, p2, y_bar, y, v, z;
  int m = m_bar + w;
  p.SetLength(m);
  p1.SetLength(m_bar);
  p2.SetLength(w);
  y_bar.SetLength(n);
  y.SetLength(n);
  v.SetLength(n);
  z.SetLength(w);
  x.SetLength(m);
  for (int i = 0; i < m; i++) {
    p[i] = (int) gsl_ran_gaussian(r, pow(log(n),0.5));
    if (i < m_bar) {
      p1[i] = p[i];
    } else {
      p2[i-m_bar] = p[i];
    }
  }
  y_bar = A_bar * (p1 - R*p2);
  y = G * p2;
  v = u - y_bar - y;
  //cout << "v: ========== " << endl;
  //cout << v << endl;
  SamplePreImageG(z, G, g, n, k, v);
  //cout << "G*z: ========== " << endl;
  //cout << G*z << endl;
  x = p + RI * z;
}

void GenerateKeys()
{
	//generate a and t
  // Primitive matrix G
  mat_ZZ G, A_bar, H, R, RI, At, temp_mat, t_G;
  vec_ZZ g, x;
  g.SetLength(k);
  for (int i= 0; i < k; i++) {
    g[i] = power2_ZZ(i);
  }

  long w = n*k;
  G.SetDims(n,w);
  printf("G: %ld x %ld\n", n, w);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < k; j++) {
      G[i][i*k + j] = power2_ZZ(j);
    }
  }

  long m_bar = 2*n;
  A_bar.SetDims(n,m_bar);
  printf("A bar: %ld x %ld\n", n, m_bar);
  for (int i = 0; i < n; i++) {
    A_bar[i][i] = 1;
  }
  for (int i = 0; i < n; i++) {
    for (int j = n+1; j < m_bar; j++) {
      A_bar[i][j] = RandomBnd(q);
    }
  }

  m = m_bar + w;
  R.SetDims(m_bar,w);
  RI.SetDims(m,w);
  printf("RI: %ld x %ld\n", m_bar, m);
  for (int i = 0; i < m_bar; i++) {
    for (int j = 0; j < w; j++) {
      RI[i][j] = RandomBnd(2);
      if (RI[i][j] == 0)
        RI[i][j] = -1;
      R[i][j] = RI[i][j];
    }
  }
  for (int i = 0; i < w; i++) {
    RI[m_bar+i][i] = 1;
  }

  At.SetDims(n,m);
  printf("calculating temp_mat\n");
  temp_mat.SetDims(n,w);
  int tstart = time(NULL);
	temp_mat = G - A_bar*R;
  int tend = time(NULL);
  printf("%d seconds to multiply!\n", tend - tstart);
  printf("At: %ld x %ld\n", n, m);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m_bar; j++) {
      At[i][j] = A_bar[i][j];
    }
  }
  printf("At, part 2\n");
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < w; j++) {
      At[i][m_bar + j] = temp_mat[i][j];
      At[i][m_bar + j] %= q;
    }
  }

  printf("Sampling\n");
  // sample x s.t. <g,x> = 0 mod q
  // g = (1, 2, 4, ... 2^{k-1})
  x.SetLength(k);

  a.SetDims(m, n);
  a = transpose(At);

  t.SetDims(m, m);
  ZZ prod_indv;
  mat_ZZ prod;
  for (int i = 0; i < m; i++) {
    vec_ZZ x;
    vec_ZZ zero;
    zero.SetLength(n);
    SamplePreImageA(x, G, A_bar, R, RI, g, n, k, m_bar, w, zero);
    for (int j = 0; j < m; j++) {
      t[j][i] = x[j];
    }
    ZZ det;
    mat_ZZ temp;
    temp = t;
    cout << "rank: " << image(det, temp) << endl;
    cout << temp << endl;
    //cout << det << endl;
    //cout << endl;
    //cout << At*x << endl;
  }
  t = transpose(t);
  prod = t*a;
  //cout << prod << endl;
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < n; j++) {
      prod[i][j] %= q;
    }
  }
  //cout << prod << endl;
  /*prod = G*t_G;
  t = transpose(RI*t_G);
  prod = t*a;
  printf("matrix is 0? %s\n", IsZero(prod) ? "yes" : "no");
  cout << "===== G" << endl << G << endl;
  cout << "===== t_G" << endl << t_G << endl;
  cout << "===== A_bar" << endl << A_bar << endl;
  cout << "===== RI" << endl << RI << endl;
  cout << "===== A" << endl << a << endl;
  cout << "===== t" << endl << t << endl;
  cout << "===== prod" << endl << prod << endl;*/
  cout << "===== t" << endl << t << endl;
  ZZ det;
  mat_ZZ t_inv;
  t_inv.SetDims(m, m);
  inv(det, t_inv, t);
  cout << "===== t_inv" << endl << t_inv << endl;
  cout << "===== t*t_inv" << endl << t*t_inv << endl;
  cout << "===== det" << endl << det << endl;

  /*ZZ_p::init(q);
  mat_ZZ_p t_p;
  t_p.SetDims(m, m);
  for (int i = 0; i < m; i++) {
    for (int j = 0; j < m; j++) {
      ZZ_p temp;
      conv(temp, t[i][j]);
      t_p[i][j] = temp;
    }
  }*/
  //mat_ZZ_p t_inv_p;
  //cout << "===== t_p" << endl << t_p << endl;
  //t_inv_p = inv(t_p);
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
