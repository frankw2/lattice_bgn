#include "bgn.h"

// prime represents the prime we are going to mod by for ciphertext
ZZ prime;

//number of total records in the database (how big the database is)
ZZ num_Records; 

//the n,m parameter in the paper
long n,m;

//the public and private keys
mat_ZZ a,t;

void GenerateKeys()
{
	//generate a and t
	
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