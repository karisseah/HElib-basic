//
// Created by karis on 6/12/18.
//

#include "encoding.h"
#include <iostream>
#include <FHE.h>

using namespace std;
using namespace NTL;

vector<vector<double>> Frac_Part(vector<vector<double>> mat) {

    int matrow = mat.size();
    int matcol = mat[0].size();

    // The remaining decimal value: initial value - integer.
    for (int i = 0; i < matrow; i++) {
        for (int j = 0; j < matcol; j++) {
            mat[i][j] = mat[i][j] - trunc(mat[i][j]);
        }
    }

    cout << "This is the fractional value: ";
    return mat;

}

vector<vector<double>> Int_Part (vector<vector<double>> mat) {

    vector<vector<double>> inv_int;

    inv_int = mat;

    int matrow = mat.size();
    int matcol = mat[0].size();

    for (int i = 0; i < matrow; i++) {
        for (int j = 0; j < matcol; j++) {
            inv_int[i][j] = trunc(mat[i][j]);
            if (inv_int[i][j] == -0) {
                inv_int[i][j] = 0;
            }
        }
    }

    // Matrix rounded down to integers.
    cout << "This is the integer value: ";
    return inv_int;

}

// integer to binary for SINGLE value.
ZZX Encode(int z) {

    ZZX ptxt;

    int remainder, digits = 0, dividend = abs(z);

    // Stops when the dividend becomes 0.
    while (dividend != 0) {
        dividend = dividend / 2;
        digits++;
    }

    int arr[digits];

    // Initialize dividend to be z again so that it doesnt use the updated dividend values.
    dividend = abs(z);

    // Array placement starts from 0. First placement is the constant's binary.
    for (int i = 0; i < digits; i++) {
        remainder = dividend % 2;
        arr[i] = remainder;
        dividend = dividend / 2;
    }

    // Prints out the array of binary numbers.
    for (int i = 0; i < digits; i++) {
        if (arr[i] == 1) {
            if (z >= 0) {
                SetCoeff(ptxt, i);
            }
            else {
                SetCoeff(ptxt, i, -1);
            }
        }
    }

    return ptxt;
}

// frac to binary for SINGLE value.
ZZX frac_encoder(double z, int phim) {

    ZZX ptxt;
    vector<double> temp;
    double frac_pt;
    double temp_z = abs(z);

    while (frac_pt != 0) {
        double store = temp_z * 2;
        double int_pt = floor(store);
        temp.push_back(int_pt);
        frac_pt = store - int_pt;
        temp_z = frac_pt;
    }

    // Add n to each exponent and flip the sign of each terms.
    for (int i = 0; i < phim; i++) {
        if (temp[i] == 1) {
            if (z >= 0) {
                SetCoeff(ptxt, (-i - 1) + phim, -1);
            }
            else {
                SetCoeff(ptxt, (-i - 1) + phim, 1);
            }
            // If fractional part overflows into int part, equate to 0.
            for (int k = 0; k <= 1*phim/7; k++) {
                ptxt[k] = 0;
            }
        }
    }

    return ptxt;

}

// Adding the int part tgt with the frac part.
vector<vector<ZZX>> frac_to_ZZX(vector<vector<double>> mat, vector<vector<double>> dec, int phim) {

    ZZX msg1;
    ZZX msg2;
    ZZX int_and_frac;
    vector<vector<ZZX>> polyn;

    int matrow = mat.size();
    int matcol = mat[0].size();

    //cout << "fraction to ZZX: " << endl;
    for (int i = 0; i < matrow; i++) {
        vector<ZZX> temp;
        for (int j = 0; j < matcol; j++) {
            double z1 = mat[i][j];
            // Int part --> ZZX.
            msg1 = Encode(z1);
            double z2 = dec[i][j];
            // Fractional part --> ZZX.
            msg2 = frac_encoder(z2, phim);
            // Adding the int and fractional parts tgt.
            int_and_frac = msg1 + msg2;
            temp.push_back(int_and_frac);
        }
        polyn.push_back(temp);
    }

    return polyn;

}

vector<vector<int>> Decrypt(FHESecKey& secretKey, vector<vector<Ctxt>> ctxt_mat, int phim) {

    ZZX temp_store_zzx;
    vector<vector<int>> vvint;

    int ctxtmat_row = ctxt_mat.size();
    int ctxtmat_col = ctxt_mat[0].size();

    for (int i = 0; i < ctxtmat_row; i++) {
        for (int j = 0; j < ctxtmat_col; j++) {
            // Decrypt each polynomial and store it in temp_store_zzx.
            secretKey.Decrypt(temp_store_zzx, ctxt_mat[i][j]);
            vector<int> vec_int;
            int integer;
            for (int k = 0; k < phim; k++) {
                // Convert each zz in temp_store_zzx into integer.
                conv(integer, temp_store_zzx[k]);
                vec_int.push_back(integer);
            }
            vvint.push_back(vec_int);
        }
    }

    return vvint;

}

vector<vector<double>> Decode(vector<vector<int>> matrix, int phim, long p) {

    int matrow = matrix.size();
    int matcol = matrix[0].size();

    vector<double> result_vec;
    vector<vector<double>> result_vvec;

    double result1 = 0;
    double result2 = 0;

    for (int i = 0; i < matrow; i++) {
        double result = 0;
        for (int j = 0; j < matcol; j++) {
            if (matrix[i][j] > floor(p/2)) {
                matrix[i][j] = matrix[i][j] - p;
            }
            // Fraction part.
            if (j > 1*phim/7) {
                result1 += (-matrix[i][j]) * pow(2, j-phim);
            }
            // Integer part.
            else {
                result2 += matrix[i][j] * pow(2, j);
            }
        }
        result = result1 + result2;
        result = fmod(result, p);
        result_vec.push_back(result);
        result1 = 0;
        result2 = 0;
    }
    result_vvec.push_back(result_vec);

    cout << "Plaintext: ";

    return result_vvec;

}