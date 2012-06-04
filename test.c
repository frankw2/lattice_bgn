#include "bgn.h"
#include <string>

void RandMat(mat_ZZ & mat, ZZ upper) {
  mat[0][0] = RandomBnd(upper);
  /*for (int i = 0; i < mat.NumRows(); i++) {
    for (int j = 0; j < mat.NumCols(); j++) {
      mat[i][j] = RandomBnd(upper);
    }
  }*/
}

void Check(mat_ZZ & one, mat_ZZ & two, string name) {
  if (one == two) {
    cout << "YAY!!!" << endl;
  } else {
    cout << name << " DOESN'T MATCH ======================" << endl;
    cout << "orig: " << endl << one << endl;
    cout << "decr: " << endl << two << endl;
  }
}

// encrypt
void EncryptAll(mat_ZZ enc_arr[], ZZ arr[], int num) {
  mat_ZZ pt;
  pt.SetDims(m, 1);
  for (int i = 0; i < num; i++) {
    pt[0][0] = arr[i];
    Encrypt(enc_arr[i], pt);
  }
}

// calculate sum
void CalculateSum(mat_ZZ & sum, mat_ZZ enc_arr[], int num) {
  if (num == 0) {
    return;
  }
  sum.SetDims(m, 1);
  for (int i = 0; i < num; i++) {
    Sum(sum, sum, enc_arr[i]);
  }
}

// calculate variance
void CalculateVar(mat_ZZ & var, mat_ZZ enc_arr[], int num) {
  if (num == 0) {
    return;
  }
  mat_ZZ sq_of_sum;
  sq_of_sum.SetDims(m, 1);
  mat_ZZ sum_of_sq;
  sum_of_sq.SetDims(m, m);
  mat_ZZ temp;
  for (int i = 0; i < num; i++) {
    Sum(sq_of_sum, sq_of_sum, enc_arr[i]);
    Multiply(temp, enc_arr[i], enc_arr[i]);
    Sum(sum_of_sq, sum_of_sq, temp);
  }
  Multiply(sq_of_sum, sq_of_sum, sq_of_sum);

  mat_ZZ n_sum_of_sq;
  n_sum_of_sq = num * sum_of_sq;
  for (int i = 0; i < n_sum_of_sq.NumRows(); i++) {
    for (int j = 0; j < n_sum_of_sq.NumCols(); j++) {
      n_sum_of_sq[i][j] %= q;
    }
  }
  mat_ZZ dec_n_sum_of_sq;
  mat_ZZ dec_sum_of_sq;
  mat_ZZ dec_sq_of_sum;
  cout << "- decrypting sum_of_sq: ";
  //cout << sum_of_sq << endl;
  DecryptProd(dec_sum_of_sq, sum_of_sq);
  cout << dec_sum_of_sq[0][0] << endl;
  cout << "- decrypting num* sum_of_sq: ";
  //cout << n_sum_of_sq << endl;
  DecryptProd(dec_n_sum_of_sq, n_sum_of_sq);
  cout << dec_n_sum_of_sq[0][0] << endl;
  cout << "- decrypting sq_of_sum: ";
  DecryptProd(dec_sq_of_sum, sq_of_sum);
  cout << dec_sq_of_sum[0][0] << endl;
  var = num*sum_of_sq - sq_of_sum;
}

int main()
{
  Init();
  GenerateKeys();

  int num = 100;
  ZZ upper = to_ZZ(sqrt(to_RR(num_Records/num)))/num;
  ZZ arr[num];
  for (int i = 0; i < num; i++) {
    arr[i] = RandomBnd(upper);
  }

  mat_ZZ enc_arr[num];
  // encrcypt
  cout << "Encrypting" << endl;
  EncryptAll(enc_arr, arr, num);
  
  // calculate sum
  cout << "Sum" << endl;
  ZZ sum;
  for (int i = 0; i < num; i++) {
   sum += arr[i];
  }
  cout << "real sum = " << sum << endl;
  mat_ZZ enc_sum;
  mat_ZZ dec_sum;
  cout << "- calculating sum" << endl;
  CalculateSum(enc_sum, enc_arr, num);
  cout << "- decrypting sum" << endl;
  Decrypt(dec_sum, enc_sum);
  cout << "encr sum = " << dec_sum[0][0] << endl;

  // calculate variance
  cout << "num^2 * Var" << endl;
  ZZ sum_of_sq, sq_of_sum, var;
  for (int i = 0; i < num; i++) {
    sq_of_sum += arr[i];
    sum_of_sq += arr[i]*arr[i];
  }
  sq_of_sum *= sq_of_sum;
  var = num*sum_of_sq - sq_of_sum;
  cout << "real num^2 * var = num * " << sum_of_sq << " - " << sq_of_sum << " = " << var << endl;
  mat_ZZ enc_var;
  mat_ZZ dec_var;
  cout << "- calculating var" << endl;
  CalculateVar(enc_var, enc_arr, num);
  cout << "- decrypting var" << endl;
  DecryptProd(dec_var, enc_var);
  cout << "encr var = " << dec_var[0][0] << endl;

  mat_ZZ pt1, pt2, pt_sum, pt_prod;
  mat_ZZ ct1, ct2, ct_sum, ct_prod;
  mat_ZZ pt1_check, pt2_check, pt_sum_check, pt_prod_check;

  pt1.SetDims(m, 1);
  pt2.SetDims(m, 1);
  pt_sum.SetDims(m, 1);
  pt_prod.SetDims(m, m);
  ct1.SetDims(m, 1);
  ct2.SetDims(m, 1);
  ct_sum.SetDims(m, 1);
  ct_prod.SetDims(m, m);
  pt1_check.SetDims(m, 1);
  pt2_check.SetDims(m, 1);
  pt_sum_check.SetDims(m, 1);
  pt_prod_check.SetDims(m, m);

  RandMat(pt1, to_ZZ(sqrt(to_RR(num_Records))));
  RandMat(pt2, to_ZZ(sqrt(to_RR(num_Records))));
  pt_sum = pt1 + pt2;
  pt_prod = pt1 * transpose(pt2);

  cout << "Encrypting" << endl;
  Encrypt(ct1, pt1);
  Encrypt(ct2, pt2);
  cout << "Summing" << endl;
  Sum(ct_sum, pt1, pt2);
  cout << "Multiplying" << endl;
  Multiply(ct_prod, pt1, pt2);
  cout << "Decrypting pt1" << endl;
  Decrypt(pt1_check, ct1);
  cout << "Decrypting pt2" << endl;
  Decrypt(pt2_check, ct2);
  cout << "Decrypting sum" << endl;
  Decrypt(pt_sum_check, ct_sum);
  cout << "Decrypting prod" << endl;
  DecryptProd(pt_prod_check, ct_prod);
  cout << "Checking" << endl;
  Check(pt1, pt1_check, "pt1");
  Check(pt2, pt2_check, "pt2");
  Check(pt_sum, pt_sum_check, "sum");
  Check(pt_prod, pt_prod_check, "sum");

	return 1;
}
