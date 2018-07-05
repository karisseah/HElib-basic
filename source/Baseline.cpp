//
// Created by karis on 14/05/18.
//

#include <gtest/gtest.h>
#include <iostream>
#include <random>
#include <FHE.h>
#include <EncryptedArray.h>
#include <iomanip>
#include "powerful.h"
#include "encoding.h"
#include "MatrixUtility.h"
#include <math.h>

using namespace std;
using namespace NTL;
using namespace chrono;

typedef high_resolution_clock Clock;

int main() {
    cout << "Hello, World!" << endl;

    cout << "-------------------- Initialization --------------------" << endl;
    auto begin_init = Clock::now();

    long m = 512;                                         // Specific modulus
    long p = 17;                                        // Plaintext base
    long r = 1;                                         // Lifting
    long L = 10;                                         // Number of levels in the modulus chain
    long c = 2;                                         // Number of columns in key-switching matrix
    long w = 2;                                         // Hamming weight of secret key
    long d = 1;                                         // Degree of the field extension
    long k = 80;                                        // Security parameter
    long s = 2;                                         // Minimum number of slots

    // Finding m
//    m = FindM(k, L, c, p, d, s, 0);                     // Find a value for m given the params
    cout << "m = " << m << endl;

    // Initializing context
    FHEcontext context(m, p, r);                        // Initialize context
    buildModChain(context, L, c);                       // Modify the context, adding primes to the modulus chain

    // Creating polynomial
    ZZX G = context.alMod.getFactorsOverZZ()[0];       // Creates the polynomial for a plaintext slot
    cout << "degree: " << deg(G) << endl;
    cout << "G: " << G << endl;

    // Generating keys
    FHESecKey secretKey(context);                       // Construct a secret key structure
    FHEPubKey &publicKey = secretKey;                   // An "upcast": FHESecKey is a subclass of FHEPubKey - Creates publicKey from secretKey
    secretKey.GenSecKey(w);                             // Generate a secret key with Hamming weight w
    addSome1DMatrices(secretKey);                       // apparently for key switching
//	addFrbMatrices(secretKey);                          // for Ctxt rotate

    // Helper Class for encryption and decryption
    EncryptedArray ea(context, G);                      // Construct EncryptedArray object ea, associated with context and G
    long nslots = ea.size();                            // Number of slots in the plaintext encoding
    cout << "slots: " << nslots << endl;

    auto end_init = Clock::now();
    cout << "FHE Ready!" << endl;
    cout << "It took: " << duration_cast<seconds>(end_init - begin_init).count() << " seconds." << endl;

    cout << endl;

    int x, y;

    ifstream myfile;
    myfile.open("/home/karis/CLionProjects/HElib-basic/dec.txt");

    // Tests if the file opens successfully.
    if (!myfile.is_open()) {
        cout << "File failed to open!" << endl;
    }

    myfile >> x;
    myfile >> y;

    vector<vector<double>> mat1;

    mat1.resize(x);
    for (int i = 0; i < mat1.size(); i++) {
        mat1[i].resize(y);
    }

    vector<vector<double>> mat2;

    mat2.resize(x);
    for (int i = 0; i < mat2.size(); i++) {
        mat2[i].resize(1);
    }

    for (int i = 0; i < x; i++) {
        mat1[i][0] = 1;
        for (int j = 1; j < y; j++) {
            myfile >> mat1[i][j];
        }
        myfile >> mat2[i][0];
    }

    cout << "Matrix X: " << mat1 << endl;
    cout << "Matrix Y: " << mat2 << endl;

    // Finding X transpose.
    vector<vector<double>> mat_trans;
    mat_trans = matrix_transpose(mat1);
    cout << "X transpose: " << mat_trans << endl;

    // Product of X transpose and X.
    vector<vector<double>> product = dotprod(mat_trans, mat1);
    cout << "Xtrans_X: " << product << endl;

    // XtransX_Inv.
    vector<vector<double>> XtransX_Inv = Inv(product);
    cout << "Inverse matrix of Xtrans_X is: " << XtransX_Inv << endl;

    // XtransX_Inv * Xtrans.
    vector<vector<double>> XtransXInv_Xtrans = dotprod(XtransX_Inv, mat_trans);

    // XtransX_Inv * Xtrans * Y.
    vector<vector<double>> matrix = dotprod(XtransXInv_Xtrans, mat2);
    cout << "This is the product matrix: " << matrix << '\n' << endl;

//    cout << "-------------------- Encryption --------------------" << endl;
//    auto begin_encrypt = Clock::now();

    // To get phim.
    int phim = to_int(context.zMStar.getPhiM());
    //cout << phim << '\n' << endl;

    // FOR XTRANSX_INV.

    cout << "FOR XTRANSX_INV: " << endl;
    cout << "This is the decimal matrix: " << XtransX_Inv << endl;

    // Fractional part.
    vector<vector<double>> inv_dec = Frac_Part(XtransX_Inv);
    cout << inv_dec << endl;

    // Integer part.
    vector<vector<double>> inv_int = Int_Part (XtransX_Inv);
    cout << inv_int << endl;

    // Encode: Conversion of fractions (int part + dec part) to ZZX.
    vector<vector<ZZX>> inv_matrix = frac_to_ZZX(inv_int, inv_dec, phim);
    cout << "encode: " << inv_matrix << '\n' << endl;

//    // Encrypt.
//    vector<vector<Ctxt>> inv_encrypt = Encrypt(publicKey, inv_matrix);
//    cout << "encrypt: " << inv_encrypt << endl;
//
//    // Decrypt.
//    vector<vector<int>> inv_decrypt = Decrypt(secretKey, inv_encrypt, phim);
//    cout << "decrypt: " << inv_decrypt << endl;
//
//    // Decode.
//    vector<vector<double>> inv_decode = Decode(inv_decrypt, phim, p);
//    cout << "decode: " << inv_decode << '\n' << endl;

    // FOR Xtrans.

    cout << "FOR X TRANSPOSE: " << endl;
    cout << "This is the decimal matrix: " << mat_trans << endl;

    // Fractional part.
    vector<vector<double>> x_dec = Frac_Part(mat_trans);
    cout << x_dec << endl;

    // Integer part.
    vector<vector<double>> x_int = Int_Part(mat_trans);
    cout << x_int << endl;

    // Encode: Conversion of fractions (int part + dec part) to ZZX.
    vector<vector<ZZX>> x_matrix = frac_to_ZZX(x_int, x_dec, phim);
    cout << "encode: " << x_matrix << '\n' << endl;

/*
    // Decode.
    vector<vector<double>> x_decode = Decode(y, x, x_enc_decrypt, phim);
    cout << "decode: " << x_decode << '\n' << endl;
*/
    // FOR Y.

    cout << "FOR Y: " << endl;
    cout << "This is the decimal matrix: " << mat2 << endl;

    // Fractional part.
    vector<vector<double>> y_dec = Frac_Part(mat2);
    cout << y_dec << endl;

    // Integer part.
    vector<vector<double>> y_int = Int_Part(mat2);
    cout << y_int << endl;

    // Encode: Conversion of fractions (int part + dec part) to ZZX.
    vector<vector<ZZX>> y_matrix = frac_to_ZZX(y_int, y_dec, phim);
    cout << "encode: " << y_matrix << '\n' << endl;

/*
    // Decode.
    vector<vector<double>> y_decode = Decode(y_enc_decrypt, phim);
    cout << "decode: " << y_decode << '\n' << endl;
*/
//    // Multiplication of enc(XtransX_Inv), enc(Xtrans) and enc(Y).
//    vector<vector<Ctxt>> inv_enc = Encrypt(publicKey, inv_matrix);
////    cout << inv_enc << '\n' << endl;
//    vector<vector<Ctxt>> xtrans_enc = Encrypt(publicKey, x_matrix);
////    cout << x_enc << '\n' << endl;
//    vector<vector<Ctxt>> y_enc = Encrypt(publicKey, y_matrix);
////    cout << y_enc << '\n' << endl;


//    cout << inv_enc.size() << endl;
//    cout << inv_enc[0].size() << endl;
////
//    cout << xtrans_enc.size() << endl;
//    cout << xtrans_enc[0].size() << endl;
//
//    cout << y_enc.size() << endl;
//    cout << y_enc[0].size() << endl;

    cout << "-------------------- Encryption --------------------" << endl;
    auto begin_encrypt = Clock::now();


    vector<vector<Ctxt>> inv_enc;
    for (int i = 0; i < inv_matrix.size(); i++) {
        vector<Ctxt> temp_ctxt;
        for (int j = 0; j < inv_matrix[0].size(); j++) {
            Ctxt enc(publicKey);
            publicKey.Encrypt(enc, inv_matrix[i][j]);
            temp_ctxt.push_back(enc);
        }
        inv_enc.push_back(temp_ctxt);
    }

    vector<vector<Ctxt>> xtrans_enc;
    for (int i = 0; i < x_matrix.size(); i++) {
        vector<Ctxt> temp_ctxt;
        for (int j = 0; j < x_matrix[0].size(); j++) {
            Ctxt enc(publicKey);
            publicKey.Encrypt(enc, x_matrix[i][j]);
            temp_ctxt.push_back(enc);
        }
        xtrans_enc.push_back(temp_ctxt);
    }

    vector<vector<Ctxt>> y_enc;
    for (int i = 0; i < y_matrix.size(); i++) {
        vector<Ctxt> temp_ctxt;
        for (int j = 0; j < y_matrix[0].size(); j++) {
            Ctxt enc(publicKey);
            publicKey.Encrypt(enc, y_matrix[i][j]);
            temp_ctxt.push_back(enc);
        }
        y_enc.push_back(temp_ctxt);
    }

    auto end_encrypt = Clock::now();
    cout << "Encryption Over!" << endl;
    cout << "It took: " << duration_cast<seconds>(end_encrypt - begin_encrypt).count() << " seconds." << endl;

    cout << "-------------------- Operation --------------------" << '\n' << endl;

    // 1st multiplication.
    vector<vector<Ctxt>> ctxt_mat_temp = mat_mat_mult(inv_enc, xtrans_enc, publicKey, secretKey);

    // 2nd multiplication.
    vector<vector<Ctxt>> ctxt_mat = mat_mat_mult(ctxt_mat_temp, y_enc, publicKey, secretKey);

    cout << "-------------------- Decryption --------------------" << endl;

    // Decrypt.
    cout << ctxt_mat.size() << " " << ctxt_mat[0].size() << endl;
    vector<vector<int>> decrypt_mat = Decrypt(secretKey, ctxt_mat, phim);
    cout << decrypt_mat << endl;

    // Decode.
    cout << Decode(decrypt_mat, phim, p) << endl;

    myfile.close();

    return 0;

}





















    //SetCoeff(msg,1,1);
    //SetCoeff(msg,2,1); // 4 = x^2

    //cout << msg << endl;

/*
    Ctxt enc(publicKey), enc2(publicKey);
    publicKey.Encrypt(enc, msg);
    //publicKey.Encrypt(enc2, msg2);

    auto end_encrypt = Clock::now();
    cout << "Encryption Over!" << endl;
    cout << "It took: " << duration_cast<seconds>(end_encrypt - begin_encrypt).count() << " seconds." << endl;


    cout << "-------------------- Operation --------------------" << endl;
    cout << "Ciphertext before operations:" << endl;
    cout << enc << endl;

    cout << "Ciphertext after addition:" << endl;
    enc.addCtxt(enc); // expect to get 4
    cout << enc << endl;

//    cout << "Ciphertext after multiplication:" << endl;
//    enc.multiplyBy(enc); // expect to get 16
//    cout << enc << endl;

    cout << "-------------------- Decryption --------------------" << endl;
    ZZX ans;
    secretKey.Decrypt(ans, enc);
    cout << "Plaintext:  " << ans << endl;
*/

