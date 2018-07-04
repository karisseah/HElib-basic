//
// Created by karis on 6/12/18.
//

#ifndef HELIB_BASIC_ENCODING_H
#define HELIB_BASIC_ENCODING_H

#include <iostream>
#include <FHE.h>

using namespace std;
using namespace NTL;

vector<vector<double>> Frac_Part(vector<vector<double>> mat);
vector<vector<double>> Int_Part(vector<vector<double>> mat);
ZZX Encode(int x);
ZZX frac_encoder(double z, int phim);
vector<vector<ZZX>> frac_to_ZZX(vector<vector<double>> mat, vector<vector<double>> dec, int phim);
vector<vector<Ctxt>> Encrypt(FHEPubKey publicKey, vector<vector<ZZX>> matrix);
vector<vector<int>> Decrypt(FHESecKey& secretKey, vector<vector<Ctxt>> ctxt_mat, int phim);
vector<vector<double>> Decode(vector<vector<int>> matrix, int phim, long p);

#endif //HELIB_BASIC_ENCODING_H