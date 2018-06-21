//
// Created by karis on 6/12/18.
//

#ifndef HELIB_BASIC_ENCODING_H
#define HELIB_BASIC_ENCODING_H

#include <iostream>
#include <FHE.h>

using namespace std;
using namespace NTL;

ZZX encode(int x);
vector<vector<double>> matrix_transpose(vector<vector<double>> product);
vector<vector<double>> dotprod(vector<vector<double>> mat1, vector<vector<double>> mat2, int x, int y);
vector<vector<double>> Inv(int y, vector<vector<double>> product);
vector<vector<ZZX>> Encrypt(long m, long p, long r, long L, long c, long w, int rows, int cols, vector<vector<double>> mat, vector<vector<double>> dec, int phim);
ZZX frac_encoder(double z, int cols, int phim);
vector<vector<ZZX>> frac_to_ZZX(int rows, int cols, vector<vector<double>> mat, vector<vector<double>> dec, int phim);
double decode(ZZX temp_store_zzx);

//vector<vector<double>> Y();
//vector<vector<double>> Xtrans_X();
//vector<vector<double>> Det(vector<vector<double>> product, int y);
//vector<vector<ZZX>> int_to_ZZX(int x, int v, vector<vector<double>> product);
//vector<vector<ZZX>> frac_to_binary (int rows, int cols, vector<vector<double>> dec, int phim);

#endif //HELIB_BASIC_ENCODING_H