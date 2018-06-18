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
vector<vector<double>> dotprod(vector<vector<int>> mat1, vector<vector<int>> mat2, int x, int y, int v);
vector<vector<ZZX>> int_to_ZZX(int x, int v, vector<vector<double>> product);
ZZX frac_encoder(double z, int cols);
vector<vector<ZZX>> frac_to_ZZX (int rows, int cols, vector<vector<double>> dec);

#endif //HELIB_BASIC_ENCODING_H