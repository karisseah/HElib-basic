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
vector<vector<double>> Frac_Part(vector<vector<double>> mat);
vector<vector<double>> Int_Part(vector<vector<double>> mat);
vector<vector<Ctxt>> Encrypt(FHEPubKey publicKey, long w, vector<vector<ZZX>> matrix);
vector<vector<int>> Decrypt(FHESecKey secretKey, vector<vector<Ctxt>> ctxt_mat, int phim);
vector<vector<double>> Decode(vector<vector<int>> matrix, int phim);
ZZX frac_encoder(double z, int cols, int phim);
vector<vector<ZZX>> frac_to_ZZX(vector<vector<double>> mat, vector<vector<double>> dec, int phim);

vector<vector<Ctxt>> test(FHEPubKey publicKey, vector<vector<ZZX>> matrix);

//int decodeINT(ZZX poly);
//vector<vector<int>> Encrypt_Decrypt(FHESecKey secretKey, long w, long L, long c, int mid, int mid1, int rows, int cols, vector<vector<ZZX>> matrix);
//vector<vector<int>> Encrypt_Decrypt(long m, long p, long r, long L, long c, long w, int rows, int cols, vector<vector<ZZX>> matrix, int phim);
//vector<vector<Ctxt>> Encrypt_Multiplication(FHESecKey secretKey, long w, int mid, int mid1, int rows, int cols, vector<vector<ZZX>> matrix1, vector<vector<ZZX>> matrix2, vector<vector<ZZX>> matrix3);
//vector<vector<Ctxt>> Encrypt_Multiplication (int mid, int mid1, long m, long p, long r, long L, long c, long w, int rows, int cols, vector<vector<ZZX>> matrix1, vector<vector<ZZX>> matrix2, vector<vector<ZZX>> matrix3);

#endif //HELIB_BASIC_ENCODING_H