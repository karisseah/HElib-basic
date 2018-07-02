//
// Created by karis on 6/12/18.
//

#include "encoding.h"
#include <iostream>
#include <FHE.h>
#include <math.h>

#include<sstream>
#include<string>

using namespace std;
using namespace NTL;
using namespace chrono;

typedef high_resolution_clock Clock;


//int decodeINT(ZZX poly)
//{
//
//    int result=0;
//    for (int i=0 ; i<deg(poly) ; i++)
//    {
//
//        int temp;
//        conv(temp,poly[i]);
//        result += temp * pow(2,i);
//
//        cout << result << endl;
//    }
//
//    return result;
//}

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

//vector<vector<double>> matrix_transpose(vector<vector<double>> mat1) {
//
//    vector<vector<double>> product_trans(mat1[0].size(), vector<double>(mat1.size()));
//
//    for (size_t i = 0; i < mat1.size(); ++i)
//        for (size_t j = 0; j < mat1[0].size(); ++j)
//            product_trans[j][i] = mat1[i][j];
//
//    return product_trans;
//
//}

//vector<vector<double>> dotprod(vector<vector<double>> mat1, vector<vector<double>> mat2, int x, int y, int z) {
//
//    int mat1row = mat1.size();
//    int mat1col = mat1[0].size();
//
//    vector<vector<double>> mult;
//    mult.resize(x, vector<double>(y));
//
////    cout << "Dot Product of Matrices: " << endl;
//    for (int i = 0; i < x; i++) {
//        for (int j = 0; j < y; j++) {
//            mult[i][j] = 0;
//            for (int k = 0; k < z; k++) {
//                mult[i][j] += mat1[i][k] * mat2[k][j];
//            }
//        }
//    }
//
//    return mult;
//
//}

//// Inverse function.
//vector<vector<double>> Inv(int y, vector<vector<double>> mult) {
//
//    double t, div;
//
//    mult.resize(y);
//    for (int i = 0; i < mult.size(); i++) {
//        mult[i].resize(2 * y);
//    }
//
//    for(int i = 0; i < y; i++) {
//        for(int j = y; j < 2 * y; j++) {
//            if (i == j - y) {
//                mult[i][j] = 1;
//            }
//            else {
//                mult[i][j] = 0;
//            }
//        }
//    }
//
//    for(int i = 0; i < y; i++)
//    {
//        t = mult[i][i];
//        for(int j = i; j < 2 * y; j++)
//            mult[i][j]=mult[i][j]/t;
//        for(int j = 0; j < y; j++)
//        {
//            if(i!=j)
//            {
//                t=mult[j][i];
//                for(int k = 0; k < 2 * y; k++)
//                    mult[j][k] = mult[j][k] - (t * mult[i][k]);
//            }
//            div = mult[0][0];
//        }
//    }
//
//    for(int i = 0; i < y; i++) {
//        for(int j = y; j < 2 * y; j++) {
//            mult[i][j] = mult[i][j] / div;
//        }
//        for (int j = 0; j < 2 * y; j++) {
//            // Bring forward the inverse matrix.
//            mult[i][j] = mult[i][j+y];
//            if (j+y >= 2 * y) {
//                break;
//            }
//        }
//    }
//
//    // Resize the matrix to keep only the inverse matrix.
//
//    mult.resize(y);
//    for (int i = 0; i < mult.size(); i++) {
//        mult[i].resize(y);
//    }
//
//    return mult;
//
//}

vector<vector<double>> Frac_Part(vector<vector<double>> mat) {

    vector<vector<double>> inv_dec;

    inv_dec = mat;

    int matrow = mat.size();
    int matcol = mat[0].size();

    // The remaining decimal value: initial value - integer.
    for (int i = 0; i < matrow; i++) {
        for (int j = 0; j < matcol; j++) {
            inv_dec[i][j] = inv_dec[i][j] - trunc(inv_dec[i][j]);
        }
    }

    cout << "This is the fractional value: ";
    return inv_dec;

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

//vector<vector<Ctxt>> Encrypt(long m, long p, long r, long L, long c, long w, int rows, int cols, vector<vector<ZZX>> matrix) {
vector<vector<Ctxt>> Encrypt(FHEPubKey publicKey, long w, vector<vector<ZZX>> matrix) {

    cout << "-------------------- Encryption --------------------" << endl;
    auto begin_encrypt = Clock::now();

//    FHEcontext context(m, p, r);
//    buildModChain(context, L, c);
//    FHESecKey secretKey(context);
//    FHEPubKey &publicKey = secretKey;
//    secretKey.GenSecKey(w);

    // Encrypt for all ZZX in mat


    int matrow = matrix.size();
    int matcol = matrix[0].size();


    vector<vector<Ctxt>> ctxt_mat;
    for (int i = 0; i < matrow; i++) {
        vector<Ctxt> temp_ctxt;
        for (int j = 0; j < matcol; j++) {
            Ctxt enc(publicKey);
            publicKey.Encrypt(enc, matrix[i][j]);
            temp_ctxt.push_back(enc);
        }
        ctxt_mat.push_back(temp_ctxt);
    }

    auto end_encrypt = Clock::now();
    cout << "Encryption Over!" << endl;
    cout << "It took: " << duration_cast<seconds>(end_encrypt - begin_encrypt).count() << " seconds." << '\n' << endl;

    cout << "-------------------- Operation --------------------" << endl;
    cout << "Ciphertext before operations:" << endl;
//    cout << ctxt_mat << endl;

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

    return ctxt_mat;

}

//vector<vector<Ctxt>> Encrypt_Multiplication(FHEcontext context, FHEPubKey *publicKey, int mid, int mid1, long m, long p, long r, long L, long c, long w, int rows, int cols, vector<vector<ZZX>> matrix1, vector<vector<ZZX>> matrix2, vector<vector<ZZX>> matrix3) {
vector<vector<Ctxt>> Encrypt_Multiplication(FHESecKey secretKey, long w, int mid, int mid1, int rows, int cols, vector<vector<ZZX>> matrix1, vector<vector<ZZX>> matrix2, vector<vector<ZZX>> matrix3) {

    cout << "-------------------- Encryption --------------------" << endl;
    auto begin_encrypt = Clock::now();

//    FHEcontext context(m, p, r);
//    buildModChain(context, L, c);
//    FHESecKey secretKey(context);
    FHEPubKey &publicKey = secretKey;
    secretKey.GenSecKey(w);

    // Encrypt for all ZZX in mat
    Ctxt enc(publicKey);
    vector<vector<Ctxt>> ctxt1_mat;
    vector<vector<Ctxt>> ctxt2_mat;
    vector<vector<Ctxt>> ctxt3_mat;

    for (int i = 0; i < rows; i++) {
        vector<Ctxt> temp1_ctxt;
        vector<Ctxt> temp2_ctxt;
        vector<Ctxt> temp3_ctxt;
        for (int j = 0; j < cols; j++) {
            publicKey.Encrypt(enc, matrix1[i][j]);
            temp1_ctxt.push_back(enc);
//            publicKey.Encrypt(enc, matrix2[i][j]);
//            temp2_ctxt.push_back(enc);
//            publicKey.Encrypt(enc, matrix3[i][j]);
//            temp3_ctxt.push_back(enc);
            cout << "j = " << j << endl;
        }
        ctxt1_mat.push_back(temp1_ctxt);
        ctxt2_mat.push_back(temp2_ctxt);
//        ctxt3_mat.push_back(temp3_ctxt);
        cout << "i = " << i << endl;
    }

    auto end_encrypt = Clock::now();
    cout << "Encryption Over!" << endl;
    cout << "It took: " << duration_cast<seconds>(end_encrypt - begin_encrypt).count() << " seconds." << '\n' << endl;

    cout << "-------------------- Operation --------------------" << endl;
    cout << "Ciphertext before operations:" << endl;
    cout << "first ctxt: " << ctxt1_mat << endl;
    cout << "second ctxt: " << ctxt2_mat << endl;
    cout << "third ctxt: " << ctxt3_mat << endl;


    cout << "Ciphertext after multiplication:" << endl;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            for (int k = 0; k < mid; k++) {
                ctxt1_mat[i][k].multiplyBy(ctxt2_mat[k][j]);
                //ctxt1_mat[i][j].multiplyBy2(ctxt2_mat[i][j], ctxt3_mat[i][j]);
            }
        }
    }

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < 1; j++) {
            for (int k = 0; k < mid1; k++) {
                ctxt1_mat[i][k].multiplyBy(ctxt3_mat[k][j]);
                //ctxt1_mat[i][j].multiplyBy2(ctxt2_mat[i][j], ctxt3_mat[i][j]);
            }
        }
    }

    return ctxt1_mat;

}

//vector<vector<int>> Decrypt(long m, long p, long r, long L, long c, long w, int rows, int cols, vector<vector<Ctxt>> ctxt_mat, int phim) {
vector<vector<int>> Decrypt(FHESecKey secretKey, vector<vector<Ctxt>> ctxt_mat, int phim) {

//    FHEcontext context(m, p, r);
//    buildModChain(context, L, c);
//    FHESecKey secretKey(context);
//    FHEPubKey &publicKey = secretKey;
//    secretKey.GenSecKey(w);

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

// Decrypt and Decode.

//    ZZX temp_store_zzx;
//    vector<vector<int>> vvint;

//    for (int i = 0; i < rows; i++) {
//        for (int j = 0; j < cols; j++) {
//            // Decrypt each polynomial and store it in temp_store_zzx.
//            secretKey.Decrypt(temp_store_zzx, ctxt_mat[i][j]);
//            cout << "zzx: " << temp_store_zzx << endl;
//            vector<int> vec_int;
//            int integer;
//            for (int k = 0; k < phim; k++) {
//                // Convert each zz in temp_store_zzx into integer.
//                conv(integer, temp_store_zzx[k]);
//                vec_int.push_back(integer);
//            }
//            vvint.push_back(vec_int);
//        }
//    }
//
//    return vvint;


}

vector<vector<double>> Decode(vector<vector<int>> matrix, int phim, long p) {

    vector<vector<vector<int>>> storage;

    int matrow = matrix.size();
    int matcol = matrix[0].size();

    storage.resize(matrow);
    for (int i = 0; i < storage.size(); i++) {
        storage[i].resize(matcol);
    }

    int index = 0;
    for (int i = 0; i < matrow; i++) {
        for (int j = 0; j < matcol; j++) {
            storage[i][j] = matrix[index++];
        }
    }

    vector<double> result_vec;
    vector<vector<double>> result_vvec;

    // Initialize sum of fractional part and sum of int part to be 0.
    double result1 = 0;
    double result2 = 0;
    for (int i = 0; i < matrow; i++) {
        for (int j = 0; j < matcol; j++) {
            for (int k = 0; k < phim; k++) {
                // Fractional part.
                if (k > phim / 2) {
                    if (storage[i][j][k] > floor(p/2)) {
                        storage[i][j][k] = storage[i][j][k] - p;
                    }
                    // -x^(n-i) --> x^i
                    result1 += (-1*storage[i][j][k]) * pow(2, k - phim);
                }
                // Integer part.
                else {
                    result2 += storage[i][j][k] * pow(2, k);
                }
            }
            double final = result1 + result2;
            result_vec.push_back(final);
            // Reset the sum to be 0 for the next iteration.
            result1 = 0;
            result2 = 0;
        }
    }

    result_vvec.resize(matrow);
    for (int i = 0; i < result_vvec.size(); i++) {
        result_vvec[i].resize(matcol);
    }

    int index1 = 0;
    for (int i = 0; i < matrow; i++) {
        for (int j = 0; j < matcol; j++) {
            result_vvec[i][j] = result_vec[index1++];
        }
    }

    cout << "-------------------- Decryption --------------------" << endl;
    cout << "Plaintext:  ";

    return result_vvec;

}

// frac to binary for SINGLE value.
ZZX frac_encoder(double z, int cols, int phim) {

    ZZX ptxt;
    vector<double> temp;
    double frac_pt;


    while (frac_pt != 0) {
        double store = z * 2;
        double int_pt = floor(store);
        temp.push_back(int_pt);
        frac_pt = store - int_pt;
        z = frac_pt;
    }

    // Add n to each exponent and flip the sign of each terms.
    for (int i = 0; i < cols + phim; i++) {
        if (temp[i] == 1) {
            SetCoeff(ptxt, (-i-1) + phim, -1);
            // If fractional part overflows into int part, equate to 0.
            for (int k = 0; k <= phim / 2; k++) {
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
        vector<ZZX> tempp;
        for (int j = 0; j < matcol; j++) {
            double z1 = mat[i][j];
            // Int part --> ZZX.
            msg1 = encode(z1);
            double z2 = dec[i][j];
            // Fractional part --> ZZX.
            msg2 = frac_encoder(z2, matcol, phim);
            // Adding the int and fractional parts tgt.
            int_and_frac = msg1 + msg2;
            tempp.push_back(int_and_frac);
        }
        polyn.push_back(tempp);
    }

    return polyn;

}














/*
// Integers to ZZX in matrix.
vector<vector<ZZX>> int_to_ZZX(int x, int v, vector<vector<double>> product) {

    ZZX msg1;

    vector<vector<ZZX>> mat_int;

    // int --> ZZX for each elements in product matrix.
    cout << "int to ZZX: " << endl;
    for (int i = 0; i < x; i++) {
        vector<ZZX> temp;
        for (int j = 0; j < v; j++) {
            int z = product[i][j];
            msg1 = encode(z);
            temp.push_back(msg1);
        }
        mat_int.push_back(temp);
    }

    return mat_int;

}

// Fractions to binary in matrix.
vector<vector<ZZX>> frac_to_binary(int rows, int cols, vector<vector<double>> dec, int phim) {


    ZZX msg2;

    vector<vector<ZZX>> binary;

    // frac --> ZZX for each elements in product matrix.
    cout << "frac to binary: " << endl;
    for (int i = 0; i < rows; i++) {
        vector<ZZX> temp;
        for (int j = 0; j < cols; j++) {
            long double z = dec[i][j];
            msg2 = frac_encoder(z, cols, phim);
            temp.push_back(msg2);
        }
        binary.push_back(temp);
    }

    return binary;

}*/

