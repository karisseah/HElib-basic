//
// Created by karis on 6/29/18.
//

#ifndef HELIB_BASIC_MATRIXUTILITY_H
#define HELIB_BASIC_MATRIXUTILITY_H

#include <iostream>
#include <FHE.h>

using namespace std;
using namespace NTL;

vector<vector<double>> matrix_transpose(vector<vector<double>> product);
vector<vector<double>> dotprod(vector<vector<double>> mat1, vector<vector<double>> mat2);
vector<vector<double>> Inv(vector<vector<double>> product);
vector<vector<Ctxt>> mat_mat_mult(vector<vector<Ctxt>> mat1, vector<vector<Ctxt>> mat2, FHEPubKey& publicKey);

#endif //HELIB_BASIC_MATRIXUTILITY_H
