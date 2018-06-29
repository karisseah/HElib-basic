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
vector<vector<Ctxt>> Encrypt(FHESecKey secretKey, long w, int rows, int cols, vector<vector<ZZX>> matrix);
vector<vector<double>> Frac_Part(vector<vector<double>> mat, int rows, int cols);
vector<vector<double>> Int_Part(vector<vector<double>> mat, int rows, int cols);
vector<vector<int>> Decrypt(FHESecKey secretKey, int rows, int cols, vector<vector<Ctxt>> ctxt_mat, int phim);
vector<vector<double>> Decode(int rows, int cols, vector<vector<int>> matrix, int phim);
ZZX frac_encoder(double z, int cols, int phim);
vector<vector<ZZX>> frac_to_ZZX(int rows, int cols, vector<vector<double>> mat, vector<vector<double>> dec, int phim);

// Shifted to MU.
//vector<vector<double>> matrix_transpose(vector<vector<double>> product);
//vector<vector<double>> dotprod(vector<vector<double>> mat1, vector<vector<double>> mat2, int x, int y, int z);
//vector<vector<double>> Inv(int y, vector<vector<double>> product);

// MISC.
vector<vector<int>> Encrypt_Decrypt(FHESecKey secretKey, long w, long L, long c, int mid, int mid1, int rows, int cols, vector<vector<ZZX>> matrix);
//vector<vector<int>> Encrypt_Decrypt(long m, long p, long r, long L, long c, long w, int rows, int cols, vector<vector<ZZX>> matrix, int phim);
vector<vector<Ctxt>> Encrypt_Multiplication(FHESecKey secretKey, long w, int mid, int mid1, int rows, int cols, vector<vector<ZZX>> matrix1, vector<vector<ZZX>> matrix2, vector<vector<ZZX>> matrix3);
//vector<vector<Ctxt>> Encrypt_Multiplication (int mid, int mid1, long m, long p, long r, long L, long c, long w, int rows, int cols, vector<vector<ZZX>> matrix1, vector<vector<ZZX>> matrix2, vector<vector<ZZX>> matrix3);

#endif //HELIB_BASIC_ENCODING_H