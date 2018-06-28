//
// Created by karis on 6/12/18.
//

#ifndef HELIB_BASIC_ENCODING_H
#define HELIB_BASIC_ENCODING_H

#include <iostream>
#include <FHE.h>

using namespace std;
using namespace NTL;

//int decodeINT(ZZX poly);
ZZX encode(int x);
vector<vector<double>> matrix_transpose(vector<vector<double>> product);
vector<vector<double>> dotprod(vector<vector<double>> mat1, vector<vector<double>> mat2, int x, int y, int z);
vector<vector<double>> Inv(int y, vector<vector<double>> product);
vector<vector<double>> Frac_Part(vector<vector<double>> mat, int rows, int cols);
vector<vector<double>> Int_Part(vector<vector<double>> mat, int rows, int cols);
//vector<vector<int>> Encrypt_Decrypt(long m, long p, long r, long L, long c, long w, int rows, int cols, vector<vector<ZZX>> matrix, int phim);
vector<vector<Ctxt>> Encrypt(long m, long p, long r, long L, long c, long w, int rows, int cols, vector<vector<ZZX>> matrix);
vector<vector<int>> Decrypt(long m, long p, long r, long L, long c, long w, int rows, int cols, vector<vector<Ctxt>> ctxt_mat, int phim);
vector<vector<double>> Decode(int rows, int cols, vector<vector<int>> matrix, int phim);
ZZX frac_encoder(double z, int cols, int phim);
vector<vector<ZZX>> frac_to_ZZX(int rows, int cols, vector<vector<double>> mat, vector<vector<double>> dec, int phim);

#endif //HELIB_BASIC_ENCODING_H