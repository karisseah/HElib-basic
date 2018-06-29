//
// Created by karis on 6/29/18.
//

#ifndef HELIB_BASIC_MATRIXUTILITY_H
#define HELIB_BASIC_MATRIXUTILITY_H

#include <iostream>
#include <FHE.h>

using namespace std;
using namespace NTL;

// Codes written but edit the rows, cols, etc.
vector<vector<double>> matrix_transpose(vector<vector<double>> product);
vector<vector<double>> dotprod(vector<vector<double>> mat1, vector<vector<double>> mat2, int x, int y, int z);
vector<vector<double>> Inv(int y, vector<vector<double>> product);
// Codes yet to write.
vector<vector<Ctxt>> mat_mat_mult();
vector<vector<Ctxt>> mat_vec_mult();

#endif //HELIB_BASIC_MATRIXUTILITY_H
