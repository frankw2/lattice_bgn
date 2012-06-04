#ifndef _bgn_h
#define _bgn_h

#include <NTL/ZZXFactoring.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/vec_ZZ.h>
#include <NTL/ZZ_pX.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_p.h>

using namespace std;
using namespace NTL;

//extern everything for testing purposes (even private values)
extern ZZ q, num_Records;
extern long n,m;
extern mat_ZZ a,t;

//generates public and private keys and also public parameters
void Init();

//generates public and private keys and also public parameters
void GenerateKeys();

//encryption function (takes in message and sets ciphertext)
void Encrypt(mat_ZZ& ciphertext, mat_ZZ& message);

//decryption function (takes in ciphertext and outputs message)
void Decrypt(mat_ZZ& message, mat_ZZ& ciphertext);

void DecryptProd(mat_ZZ& message, mat_ZZ& ciphertext);

void Sum(mat_ZZ & sum, mat_ZZ & ct1, mat_ZZ & ct2);

void Multiply(mat_ZZ & prod, mat_ZZ & ct1, mat_ZZ & ct2);

#endif
