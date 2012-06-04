#include "bgn.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <NTL/LLL.h>
#include <cmath>
#include <ctime>

// prime represents the prime we are going to mod by for ciphertext
ZZ q;

//q is 2^k
int k;

int c;
RR beta;

//number of total records in the database (how big the database is, i.e. 2^16)
ZZ num_Records; 

//the n,m parameter in the paper
long n,m;

long w, m_bar;

//the public and private keys
mat_ZZ a,t;

const gsl_rng_type * T;
gsl_rng *r;

void Init() {
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);
  gsl_rng_set(r, time(NULL));
  ZZ seed;
  seed = time(NULL);
  SetSeed(seed); 

  c = 2;
  n = 10;

  RR two = to_RR(2);
  RR cc = to_RR(c);
  RR nn = to_RR(n);
  RR q_min = power(two, 20) * power(cc+4, 3) * power(nn, 3*c+4) * power(log(nn), 5);
  double db = log(to_ZZ(q_min))/log(2);
  k = (int) db;
  q = power2_ZZ(k-1);
  q += RandomBnd(q);
  q = NextPrime(q);
  k = (int) ceil(log(q)/log(2));

  w = n*k;
  m_bar = 2*n;
  m = m_bar + w;
  RR qq = to_RR(q);
  beta =  1.0 / (27 * power(nn, 1 + 3 * c / 2) * log(nn) * log(qq) * sqrt(qq * m));

  num_Records = 1 << 16;
  cout << "n = " << n << endl;
  cout << "k = " << k << endl;
  cout << "q = " << q << endl;
  cout << "beta = " << beta << endl;
  cout << "num = " << num_Records << endl;
  printf("done with init\n");
}

#define PI 3.14159265358979
void SamplePsiFromZ(ZZ & x, RR & beta) {
  RR y = to_RR(gsl_ran_gaussian(r, to_double(beta / sqrt(2*PI))));
  x = to_ZZ(round(y * to_RR(q))) % q;
}

void SampleZ(ZZ & x, RR & s, RR & c, int k) {
  while (true) {
    ZZ lower, upper;
    RR temp;
    ceil(temp, c-s*log(k));
    lower = to_ZZ(temp);
    floor(temp, c+s*log(k));
    upper = to_ZZ(temp);
    x = RandomBnd(upper-lower+1);
    x += lower;
    //double p1 = gsl_ran_gaussian_pdf(to_double(to_RR(x)-c), to_double(s));
    RR p = exp(-PI * power(to_RR(x)-c, 2) / power(s, 2));
    RR random = random_RR();
    //cout << p << " " << random << " " << x << " " << lower << " " << upper << " " << c << " " << s << " " << s * log(k) << endl;
    if (random < p) {
      break;
    }
  }
}

void to_mat_RR(mat_RR & out, mat_ZZ in) {
  out.SetDims(in.NumRows(), in.NumCols());
  for (int i = 0; i < in.NumRows(); i++) {
    for (int j = 0; j < in.NumCols(); j++) {
      out[i][j] = to_RR(in[i][j]);
    }
  }
}

void to_vec_RR(vec_RR & out, vec_ZZ in) {
  out.SetLength(in.length());
  for (int i = 0; i < in.length(); i++) {
    out[i] = to_RR(in[i]);
  }
}

// Samples z such that G*out = target. G is n x nk.
void SamplePreImageG(vec_ZZ & out, mat_ZZ & G, mat_ZZ & Sk, vec_ZZ target) {
  mat_ZZ St;
  St = transpose(Sk);
  int ct = 0;
  for (int i = 0; i < n; i++) {
    vec_ZZ v;
    v.SetLength(w);
    v = St[ct];
    for (int j = 0; j < k; j++) {
      out[i*k+j] = v[j] + ((target[i] >> j) & 1);
    }
    ct = (ct+1)%k;
  }
}

void GenerateKeys()
{
	//generate a and t
  // Primitive matrix G
  mat_ZZ G, Sk, S, A_bar, H, R, RI, At, temp_mat, t_G;
  mat_RR Sk_tilde;
  vec_ZZ g;
  g.SetLength(k);
  for (int i= 0; i < k; i++) {
    g[i] = power2_ZZ(i);
  }

  G.SetDims(n,w);
  printf("G: %ld x %ld\n", n, w);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < k; j++) {
      G[i][i*k + j] = power2_ZZ(j);
    }
  }

  Sk.SetDims(k,k);
  printf("Sk: %d x %d\n", k, k);
  for (int i = 0; i < Sk.NumRows(); i++) {
    Sk[i][i] = 2;
    if (i > 0) {
      Sk[i][i-1] = -1;
    }
    Sk[i][Sk.NumCols()-1] = (q >> i) & 1;
  }
  //cout << Sk << endl;

  S.SetDims(w,w);
  printf("S: %ld x %ld\n", w, w);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < k; j++) {
      for (int l = 0; l < k; l++) {
        S[k*i+j][k*i+l] = Sk[j][l];
      }
    }
  }

  mat_RR mu;
  vec_RR c;
  mat_RR S_transpose;
  to_mat_RR(S_transpose, transpose(Sk));
  ComputeGS(transpose(Sk), mu, c);
  //cout << c << endl;
  Sk_tilde.SetDims(Sk.NumRows(), Sk.NumCols());
  for (int i = 0; i < Sk_tilde.NumRows(); i++) {
    vec_RR sum;
    sum.SetLength(Sk_tilde.NumCols());
    for (int j = 0; j < i; j++) {
      sum += mu[i][j] * Sk_tilde[j];
    }
    Sk_tilde[i] = S_transpose[i] - sum;
  }
  Sk_tilde = transpose(Sk_tilde);
  //cout << Sk_tilde << endl<< endl;

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

  printf("Setting A\n");
  a.SetDims(m, n);
  a = transpose(At);

  printf("Setting T\n");
  t.SetDims(m, m);
  mat_ZZ left, right, A_target, W;
  left.SetDims(m, m);
  right.SetDims(m, m);
  A_target.SetDims(n, m_bar);
  W.SetDims(w, m_bar);
  printf("- Left\n");
  // set up left
  for (int i = 0; i < m; i++) {
    left[i][i] = 1;
  }
  for (int i = 0; i < m_bar; i++) {
    for (int j = 0; j < w; j++) {
      left[i][m_bar+j] = R[i][j];
    }
  }

  printf("- Target\n");
  // set up A_target = A [I | 0]^T
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < m_bar; j++) {
      A_target[i][j] = At[i][j];
    }
  }
  A_target = transpose(A_target); // to get the columns of A_target before transpose

  printf("- Right\n");
  // set up W
  for (int i = 0; i < m_bar; i++) {
    vec_ZZ x;
    x.SetLength(w);
    vec_ZZ target;
    target = A_target[i];
    SamplePreImageG(x, G, Sk, target);
    for (int j = 0; j < w; j++) {
      W[j][i] = -x[j];
    }
  }
  // set up right
  for (int i = 0; i < m_bar; i++) {
    right[i][i] = 1;
  }
  for (int i = 0; i < w; i++) {
    for (int j = 0; j < m_bar; j++) {
      right[m_bar+i][j] = W[i][j];
    }
  }
  for (int i = 0; i < w; i++) {
    for (int j = 0; j < w; j++) {
      right[m_bar+i][m_bar+j] = S[i][j];
    }
  }
  t = left * right;
  t = transpose(t);

  cout << "Checking T * A" << endl;
  mat_ZZ prod;
  prod = t*a;
  for (int i = 0; i < prod.NumRows(); i++) {
    for (int j = 0; j < prod.NumCols(); j++) {
      prod[i][j] %= q;
    }
  }
  //cout << determinant(t) << " " << power(q, n) << endl;
  printf("determinant T is q^n? %s\n", determinant(t) == power(q, n) ? "yes" : "no");
  printf("matrix T*A is 0? %s\n", IsZero(prod) ? "yes" : "no");
}

//message should be square matrix m x m 
void Encrypt(mat_ZZ& ciphertext, mat_ZZ& message)
{
  ciphertext.SetDims(message.NumRows(), message.NumCols());
	//set and create random matrix
	mat_ZZ s;
	s.SetDims(n,message.NumCols());
	for (int i = 0; i < n; i++){
		for (int j = 0; j < ciphertext.NumCols(); j++){
			s[i][j] = RandomBnd(q);
		}
	}
  
  //cout << beta << endl;
	mat_ZZ term1;
	term1 = a*s;
	
  mat_ZZ term2;
  term2.SetDims(ciphertext.NumRows(), ciphertext.NumCols());
  for (int i = 0; i < ciphertext.NumRows(); i++) {
    for (int j = 0; j < ciphertext.NumCols(); j++) {
      SamplePsiFromZ(term2[i][j], beta);
    }
  }
  term2 = num_Records * term2;
  //cout << term2 << endl;
	ciphertext = term1 + term2 + message;
  
	for (int i = 0; i < ciphertext.NumRows(); i++){
		for (int j = 0; j < ciphertext.NumCols(); j++){
      ciphertext[i][j] %= q;
    }
  }
}

void Decrypt(mat_ZZ& message, mat_ZZ& ciphertext)
{
	mat_ZZ e;

	e = t*ciphertext;
	for (int i = 0; i < e.NumRows(); i++){
		for (int j = 0; j < e.NumCols(); j++){
			e[i][j] = e[i][j] % q;
      if (e[i][j] > q/2) {
        e[i][j] -= q;
      }
		}
	}
  ZZ det;
  mat_ZZ t_inv, t_transpose_inv;
  t_inv.SetDims(m, m);
  t_transpose_inv.SetDims(m,m);
  inv(det, t_inv, t);
  t_transpose_inv = transpose(t_inv);

	message = t_inv*e;//*t_transpose_inv;

	for (int i = 0; i < message.NumRows(); i++){
		for (int j = 0; j < message.NumCols(); j++){
      message[i][j] /= det;
      //message[i][j] /= det;
      message[i][j] %= num_Records;
		}
	}
}


void DecryptProd(mat_ZZ& message, mat_ZZ& ciphertext)
{
	mat_ZZ e;
	
	e = t*ciphertext*transpose(t);
	for (int i = 0; i < e.NumRows(); i++){
		for (int j = 0; j < e.NumCols(); j++){
			e[i][j] = e[i][j] % q;
      if (e[i][j] > q/2) {
        e[i][j] -= q;
      }
		}
	}
  //cout << e << endl;
  ZZ det;
  mat_ZZ t_inv, t_transpose_inv;
  t_inv.SetDims(m, m);
  t_transpose_inv.SetDims(m,m);
  inv(det, t_inv, t);
  t_transpose_inv = transpose(t_inv);

	message = t_inv*e*t_transpose_inv;

	for (int i = 0; i < message.NumRows(); i++){
		for (int j = 0; j < message.NumCols(); j++){
      message[i][j] /= det;
      message[i][j] /= det;
      message[i][j] %= num_Records;
		}
	}
}

void MatrixMod(mat_ZZ & mat, ZZ mod) {
  for (int i = 0; i < mat.NumRows(); i++) {
    for (int j = 0; j < mat.NumCols(); j++) {
      mat[i][j] %= mod;
    }
  }
}

void Sum(mat_ZZ & sum, mat_ZZ & ct1, mat_ZZ & ct2) {
  sum = ct1 + ct2;
  MatrixMod(sum, q);
}

void Multiply(mat_ZZ & prod, mat_ZZ & ct1, mat_ZZ & ct2) {
  prod = ct1 * transpose(ct2);
  MatrixMod(prod, q);
}
