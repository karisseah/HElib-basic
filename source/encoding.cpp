//
// Created by karis on 6/12/18.
//

#include "encoding.h"
#include <iostream>
#include <FHE.h>

using namespace std;
using namespace NTL;
using namespace chrono;

typedef high_resolution_clock Clock;

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

//    ZZX msg;

    // Prints out the array of binary numbers.
    for (int i = 0; i < digits; i++) {
        if (arr[i] == 1) {
            SetCoeff(ptxt, i);
        }
    }

    return ptxt;
}


vector<vector<int>> dotprod(vector<vector<int>> mat1, vector<vector<int>> mat2, int x, int y, int v) {

    vector<vector<int>> mult;

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

vector<vector<ZZX>> int_to_ZZX(int x, int v, vector<vector<double>> product) {

    ZZX msg;
    vector<vector<ZZX>> mat;

    // int --> ZZX for each elements in product matrix.
    cout << "int to ZZX: " << endl;
    for (int i = 0; i < x; i++) {
        vector<ZZX> temp;
        for (int j = 0; j < v; j++) {
            int z = product[i][j];
            msg = encode(z);
            temp.push_back(msg);
        }
        mat.push_back(temp);
    }

    return mat;

}

//vector<vector<ZZX>> Int_Encoder(int x, int v, FHEPubKey &publicKey, vector<vector<ZZX>> matrix, vector<vector<int>> product, FHESecKey secretKey(context)) {
//
//    auto begin_encrypt = Clock::now();
//
//    // vector<vector<int>> ---> vector< vector<ZZX>>
//    matrix = int_to_ZZX(x, v, product);
//
//    cout << matrix << endl;
//
//    // Encrypt for all ZZX in mat
//    Ctxt enc(publicKey);
//    vector<vector<Ctxt>> ctxt_mat;
//
//    for (int i = 0; i < x; i++) {
//        vector<Ctxt> temp_ctxt;
//        for (int j = 0; j < v; j++) {
//            publicKey.Encrypt(enc, matrix[i][j]);
//            temp_ctxt.push_back(enc);
//        }
//        ctxt_mat.push_back(temp_ctxt);
//    }
//
//    auto end_encrypt = Clock::now();
//    cout << "Encryption Over!" << endl;
//    cout << "It took: " << duration_cast<seconds>(end_encrypt - begin_encrypt).count() << " seconds." << '\n' << endl;
//
//    cout << "-------------------- Operation --------------------" << endl;
//    cout << "Ciphertext before operations:" << endl;
//    cout << ctxt_mat << endl;
//
////    cout << "Ciphertext after addition:" << endl;
////    for (int i = 0; i < x; i++) {
////        for (int j = 0; j < v; j++) {
////            ctxt_mat[i][j].addCtxt(ctxt_mat[i][j]);
////        }
////    }
////    cout << ctxt_mat << endl;
//
////    cout << "Ciphertext after multiplication:" << endl;
////    for (int i = 0; i < x; i++) {
////        for (int j = 0; j < v; j++) {
////            ctxt_mat[i][j].multiplyBy(ctxt_mat[i][j]);
////        }
////    }
////    cout << ctxt_mat << endl;
//
//    vector<vector<ZZX>> mat_ans;
//    ZZX temp_store;
//
//    for (int i = 0; i < x; i++) {
//        vector<ZZX> mat_temp;
//        vector<Ctxt> temp_ctxt;
//        for (int j = 0; j < v; j++) {
//            secretKey.Decrypt(temp_store, ctxt_mat[i][j]);
//            mat_temp.push_back(temp_store);
//        }
//        mat_ans.push_back(mat_temp);
//    }
//
//    cout << "-------------------- Decryption --------------------" << endl;
//    cout << "Plaintext:  " << mat_ans << endl;
//
//    return 0;
//
//}

vector<vector<ZZX>> frac_to_ZZX() {


}