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
vector<vector<double>> Y();
vector<vector<double>> Xtrans_X();
vector<vector<double>> matrix_transpose(vector<vector<double>> product);
vector<vector<double>> dotprod(vector<vector<double>> mat1, vector<vector<double>> mat2, int x, int y);
vector<vector<double>> Inv(vector<vector<double>> product, int y);
vector<vector<ZZX>> int_to_ZZX(int x, int v, vector<vector<double>> product);
vector<vector<ZZX>> Encrypt(long m, long p, long r, long L, long c, long w, int x, int v, vector<vector<double>> product);
ZZX frac_encoder(double z, int cols, int phim);
vector<vector<ZZX>> frac_to_binary (int rows, int cols, vector<vector<double>> dec, int phim);

#endif //HELIB_BASIC_ENCODING_H