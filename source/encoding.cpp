//
// Created by karis on 6/12/18.
//

#include "encoding.h"
#include <iostream>
#include <FHE.h>
#include <math.h>

using namespace std;
using namespace NTL;
using namespace chrono;

typedef high_resolution_clock Clock;

// integer to binary for SINGLE value.
ZZX encode(int z) {

    ZZX ptxt;

    int remainder, digits = 0, dividend = z;

    // Stops when the dividend becomes 0.
    while (dividend != 0) {
        dividend = dividend / 2;
        digits++;
    }

    int arr[digits];

    // Initialize dividend to be z again so that it doesnt use the updated dividend values.
    dividend = z;

    // Array placement starts from 0. First placement is the constant's binary.
    for (int i = 0; i < digits; i++) {
        remainder = dividend % 2;
        arr[i] = remainder;
        dividend = dividend / 2;
    }

    // Prints out the array of binary numbers.
    for (int i = 0; i < digits; i++) {
        if (arr[i] == 1) {
            SetCoeff(ptxt, i);
        }
    }

    return ptxt;
}

vector<vector<double>> dotprod(vector<vector<int>> mat1, vector<vector<int>> mat2, int x, int y, int v) {

    vector<vector<double>> mult;

    // Dot Product.
    mult.resize(x);
    for (int i = 0; i < mult.size(); i++) {
        mult[i].resize(v);
    }

    cout << "Dot Product of Matrices: " << endl;
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < v; j++) {
            mult[i][j] = 0;
            for (int k = 0; k < y; k++) {
                mult[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }

    return mult;

}

// Integers to ZZX in matrix.
vector<vector<ZZX>> int_to_ZZX(int x, int v, vector<vector<double>> product) {

    ZZX msg;
    vector<vector<ZZX>> mat_int;

    // int --> ZZX for each elements in product matrix.
    cout << "int to ZZX: " << endl;
    for (int i = 0; i < x; i++) {
        vector<ZZX> temp;
        for (int j = 0; j < v; j++) {
            int z = product[i][j];
            msg = encode(z);
            temp.push_back(msg);
        }
        mat_int.push_back(temp);
    }

    return mat_int;

}

vector<vector<ZZX>> Encrypt(long m, long p, long r, long L, long c, long w, int x, int v, vector<vector<double>> product) {

    FHEcontext context(m, p, r);
    buildModChain(context, L, c);
    FHESecKey secretKey(context);
    FHEPubKey &publicKey = secretKey;
    secretKey.GenSecKey(w);

    vector<vector<ZZX>> matrix;

    auto begin_encrypt = Clock::now();

    // vector<vector<int>> ---> vector< vector<ZZX>>
    matrix = int_to_ZZX(x, v, product);

    cout << matrix << endl;

    // Encrypt for all ZZX in mat
    Ctxt enc(publicKey);
    vector<vector<Ctxt>> ctxt_mat;

    for (int i = 0; i < x; i++) {
        vector<Ctxt> temp_ctxt;
        for (int j = 0; j < v; j++) {
            publicKey.Encrypt(enc, matrix[i][j]);
            // cant run from here onwards.
            temp_ctxt.push_back(enc);
        }
        ctxt_mat.push_back(temp_ctxt);
    }

    auto end_encrypt = Clock::now();
    cout << "Encryption Over!" << endl;
    cout << "It took: " << duration_cast<seconds>(end_encrypt - begin_encrypt).count() << " seconds." << '\n' << endl;

    cout << "-------------------- Operation --------------------" << endl;
    cout << "Ciphertext before operations:" << endl;
    cout << ctxt_mat << endl;

//    cout << "Ciphertext after addition:" << endl;
//    for (int i = 0; i < x; i++) {
//        for (int j = 0; j < v; j++) {
//            ctxt_mat[i][j].addCtxt(ctxt_mat[i][j]);
//        }
//    }
//    cout << ctxt_mat << endl;

//    cout << "Ciphertext after multiplication:" << endl;
//    for (int i = 0; i < x; i++) {
//        for (int j = 0; j < v; j++) {
//            ctxt_mat[i][j].multiplyBy(ctxt_mat[i][j]);
//        }
//    }
//    cout << ctxt_mat << endl;

    vector<vector<ZZX>> mat_ans;
    ZZX temp_store;

    for (int i = 0; i < x; i++) {
        vector<ZZX> mat_temp;
        vector<Ctxt> temp_ctxt;
        for (int j = 0; j < v; j++) {
            secretKey.Decrypt(temp_store, ctxt_mat[i][j]);
            mat_temp.push_back(temp_store);
        }
        mat_ans.push_back(mat_temp);
    }

    cout << "-------------------- Decryption --------------------" << endl;
//    cout << "Plaintext:  " << mat_ans << endl;
    cout << "Plaintext:  " << endl;

    return mat_ans;

}

// frac to binary for SINGLE value.
ZZX frac_encoder(double z, int cols, int phim) {

    ZZX ptxt;
    vector<double> temp;
    double frac_pt;

    while (frac_pt != 0) {
        double store = z * 2;
        int int_pt = floor(store);
        temp.push_back(int_pt);
        frac_pt = store - int_pt;
        z = frac_pt;
    }

    // Add n to each exponent and flip the sign of each terms.
    for (int i = 0; i < cols + phim; i++) {
        if (temp[i] == 1) {
            SetCoeff(ptxt, (-i-1) + phim, -1);
        }
    }

    return ptxt;

}

// Fractions to binary in matrix.
vector<vector<ZZX>> frac_to_binary(int rows, int cols, vector<vector<double>> dec, int phim) {

    ZZX msg;
    vector<vector<ZZX>> binary;

    // frac --> ZZX for each elements in product matrix.
    cout << "frac to binary: " << endl;
    for (int i = 0; i < rows; i++) {
        vector<ZZX> temp;
        for (int j = 0; j < cols; j++) {
            long double z = dec[i][j];
            msg = frac_encoder(z, cols, phim);
            temp.push_back(msg);
        }
        binary.push_back(temp);
    }

    return binary;

}
