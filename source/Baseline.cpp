//
// Created by junjie on 5/14/18.
//

#include <gtest/gtest.h>
#include <iostream>
#include <random>
#include <FHE.h>
#include <EncryptedArray.h>
#include <iomanip>
#include "powerful.h"
#include "encoding.h"
#include <math.h>

using namespace std;
using namespace NTL;
using namespace chrono;

typedef high_resolution_clock Clock;

int main() {
    cout << "Hello, World!" << endl;

    cout << "-------------------- Initialization --------------------" << endl;
    auto begin_init = Clock::now();

    long m = 190;                                         // Specific modulus
    long p = 17;                                        // Plaintext base
    long r = 1;                                         // Lifting
    long L = 5;                                         // Number of levels in the modulus chain
    long c = 1;                                         // Number of columns in key-switching matrix
    long w = 2;                                         // Hamming weight of secret key
    long d = 1;                                         // Degree of the field extension
    long k = 80;                                        // Security parameter
    long s = 0;                                         // Minimum number of slots

    // Finding m
    //m = FindM(k, L, c, p, d, s, 0);                     // Find a value for m given the params
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
//	addSome1DMatrices(secretKey);                       // apparently for key switching
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
    vector<vector<double>> product = dotprod(mat_trans, mat1, y, y, x);
    cout << "Xtrans_X: " << product << endl;

    // XtransX_Inv.
    vector<vector<double>> XtransX_Inv = Inv(y, product);
    cout << "Inverse matrix of Xtrans_X is: " << XtransX_Inv << endl;

    // XtransX_Inv * Xtrans.
    vector<vector<double>> XtransXInv_Xtrans = dotprod(XtransX_Inv, mat_trans, y, x, y);

    // XtransX_Inv * Xtrans * Y.
    vector<vector<double>> matrix = dotprod(XtransXInv_Xtrans, mat2, y, 1, x);
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
    vector<vector<double>> inv_dec = Frac_Part(XtransX_Inv, y, y);
    cout << inv_dec << endl;

    // Integer part.
    vector<vector<double>> inv_int = Int_Part (XtransX_Inv, y, y);
    cout << inv_int << endl;

    // Encode: Conversion of fractions (int part + dec part) to ZZX.
    vector<vector<ZZX>> inv_matrix = frac_to_ZZX(y, y, inv_int, inv_dec, phim);
    cout << "encode: " << inv_matrix << endl;

    // Encrypt.
    //vector<vector<Ctxt>> inv_encrypt = Encrypt(m, p, r, L, c, w, y, y, inv_matrix);
    //cout << "encrypt: " << inv_encrypt << endl;
/*
    // Encrypt and Decrypt.
    vector<vector<int>> inv_enc_decrypt = Encrypt_Decrypt(m, p, r, L, c, w, y, y, inv_matrix, phim);
    cout << "encrypt and decrypt: " << inv_enc_decrypt << endl;

    // Decrypt.
    //vector<vector<int>> inv_decrypt = Decrypt(m, p, r, L, c, w, y, y, inv_encrypt, phim);
    //cout << "decrypt: " << inv_decrypt << endl;

    // Decode.
    vector<vector<double>> inv_decode = Decode(y, y, inv_enc_decrypt, phim);
    cout << "decode: " << inv_decode << '\n' << endl;
*/
    // FOR Xtrans.

    cout << "FOR X: " << endl;
    cout << "This is the decimal matrix: " << mat_trans << endl;

    // Fractional part.
    vector<vector<double>> x_dec = Frac_Part(mat_trans, y, x);
    cout << x_dec << endl;

    // Integer part.
    vector<vector<double>> x_int = Int_Part (mat_trans, y, x);
    cout << x_int << endl;

    // Encode: Conversion of fractions (int part + dec part) to ZZX.
    vector<vector<ZZX>> x_matrix = frac_to_ZZX(y, x, x_int, x_dec, phim);
    cout << "encode: " << x_matrix << endl;

/*    // Encrypt and Decrypt.
    vector<vector<int>> x_enc_decrypt = Encrypt_Decrypt(m, p, r, L, c, w, y, x, x_matrix, phim);
    cout << "encrypt and decrypt: " << x_enc_decrypt << endl;

    // Decode.
    vector<vector<double>> x_decode = Decode(y, x, x_enc_decrypt, phim);
    cout << "decode: " << x_decode << '\n' << endl;
*/
    // FOR Y.

    cout << "FOR Y: " << endl;
    cout << "This is the decimal matrix: " << mat2 << endl;

    // Fractional part.
    vector<vector<double>> y_dec = Frac_Part(mat2, x, 1);
    cout << y_dec << endl;

    // Integer part.
    vector<vector<double>> y_int = Int_Part(mat2, x, 1);
    cout << y_int << endl;

    // Encode: Conversion of fractions (int part + dec part) to ZZX.
    vector<vector<ZZX>> y_matrix = frac_to_ZZX(x, 1, y_int, y_dec, phim);
    cout << "encode: " << y_matrix << endl;

    // Encrypt and Decrypt.
/*    vector<vector<int>> y_enc_decrypt = Encrypt_Decrypt(m, p, r, L, c, w, x, 1, y_matrix, phim);
    cout << "encrypt and decrypt: " << y_enc_decrypt << endl;

    // Decode.
    vector<vector<double>> y_decode = Decode(x, 1, y_enc_decrypt, phim);
    cout << "decode: " << y_decode << '\n' << endl;
*/
    // Multiplication of enc(XtransX_Inv), enc(Xtrans) and enc(Y).
    vector<vector<Ctxt>> inv_enc = Encrypt(m, p, r, L, c, w, y, y, inv_matrix);
    vector<vector<Ctxt>> x_enc = Encrypt(m, p, r, L, c, w, y, x, x_matrix);
    vector<vector<Ctxt>> y_enc = Encrypt(m, p, r, L, c, w, x, 1, y_matrix);




//    // FOR XTRANSX_INV.
//
//    cout << "This is the decimal matrix: " << '\n' << XtransX_Inv << '\n' << endl;
//
//    vector<vector<double>> inv_dec;
//
//    inv_dec = XtransX_Inv;
//
//    // The remaining decimal value: initial value - integer.
//    for (int i = 0; i < y; i++) {
//        for (int j = 0; j < y; j++) {
//            inv_dec[i][j] = inv_dec[i][j] - trunc(inv_dec[i][j]);
//        }
//    }
//
//    cout << "This is the fractional value: " << '\n' << inv_dec << '\n' << endl;
//
//    vector<vector<double>> inv_int;
//
//    inv_int = XtransX_Inv;
//
//    for (int i = 0; i < y; i++) {
//        for (int j = 0; j < y; j++) {
//            inv_int[i][j] = trunc(XtransX_Inv[i][j]);
//            if (inv_int[i][j] == -0) {
//                inv_int[i][j] = 0;
//            }
//        }
//    }
//
//    // Matrix rounded down to integers.
//    cout << "This is the integer value: " << '\n' << inv_int << '\n' << endl;
//
//    // Conversion of fractions (int part + dec part) to ZZX.
//    vector<vector<ZZX>> inv_matrix;
//    inv_matrix = frac_to_ZZX(y, y, inv_int, inv_dec, phim);
//
//    vector<vector<double>> Enc_inv = Encrypt_Decrypt(m, p, r, L, c, w, y, y, inv_int, inv_dec, phim, inv_matrix);
//    cout << Enc_inv << '\n' << endl;
//
//    // FOR X.
//
//    cout << "This is the decimal matrix: " << '\n' << mat1 << '\n' << endl;
//
//    vector<vector<double>> x_dec;
//
//    x_dec = mat1;
//
//    // The remaining decimal value: initial value - integer.
//    for (int i = 0; i < x; i++) {
//        for (int j = 0; j < y; j++) {
//            x_dec[i][j] = x_dec[i][j] - trunc(x_dec[i][j]);
//        }
//    }
//
//    cout << "This is the fractional value: " << '\n' << x_dec << '\n' << endl;
//
//    vector<vector<double>> x_int;
//
//    x_int = mat1;
//
//    for (int i = 0; i < x; i++) {
//        for (int j = 0; j < y; j++) {
//            x_int[i][j] = trunc(mat1[i][j]);
//        }
//    }
//
//    // Matrix rounded down to integers.
//    cout << "This is the integer value: " << '\n' << x_int << '\n' << endl;
//
//    // Conversion of fractions (int part + dec part) to ZZX.
//    vector<vector<ZZX>> x_matrix = frac_to_ZZX(x, y, x_int, x_dec, phim);
//
//    vector<vector<double>> Enc_x = Encrypt_Decrypt(m, p, r, L, c, w, x, y, x_int, x_dec, phim, x_matrix);
//    cout << Enc_x << '\n' << endl;
//
//    // FOR Y.
//
//    cout << "This is the decimal matrix: " << '\n' << mat2 << '\n' << endl;
//
//    vector<vector<double>> y_dec;
//
//    y_dec = mat2;
//
//    // The remaining decimal value: initial value - integer.
//    for (int i = 0; i < x; i++) {
//        for (int j = 0; j < 1; j++) {
//            y_dec[i][j] = y_dec[i][j] - trunc(y_dec[i][j]);
//        }
//    }
//
//    cout << "This is the fractional value: " << '\n' << y_dec << '\n' << endl;
//
//    vector<vector<double>> y_int;
//
//    y_int = mat2;
//
//    for (int i = 0; i < x; i++) {
//        for (int j = 0; j < 1; j++) {
//            y_int[i][j] = trunc(mat2[i][j]);
//        }
//    }
//
//    // Matrix rounded down to integers.
//    cout << "This is the integer value: " << '\n' << y_int << '\n' << endl;
//
//    // Conversion of fractions (int part + dec part) to ZZX.
//    vector<vector<ZZX>> y_matrix = frac_to_ZZX(x, 1, y_int, y_dec, phim);
//
//    vector<vector<double>> Enc_y = Encrypt_Decrypt(m, p, r, L, c, w, x, 1, y_int, y_dec, phim, y_matrix);
//    cout << Enc_y << '\n' << endl;




}


















/*
    // Before fractional encoder.

    int x, y, u, v;

    ifstream myfile;
    myfile.open("/home/karis/CLionProjects/HElib-basic/matrix.txt");

    // Tests if the file opens successfully.
    if (!myfile.is_open()) {
        cout << "File failed to open!" << endl;
    }

    myfile >> x;
    //cout << x << endl; // no. of rows
    myfile >> y;
    //cout << y << endl; // no. of cols

    int row1 = 0;
    vector<vector<double>> mat1;
    vector<double> rowVector1(y);

    cout << "This is the first matrix: " << endl;
    while (myfile.good()) {
        mat1.push_back(rowVector1); // add a new row,
        for (int col = 0; col < y; col++) {
            myfile >> mat1[row1][col]; // fill the row with col elements
        }
        row1++;
        if (row1 >= x) {
            break;
        }
    }
    cout << mat1 << '\n' << endl;

    myfile >> u;
    //cout << u << endl; // no. of rows
    myfile >> v;
    //cout << v << endl; // no. of cols

    assert (y == u);

    int row2 = 0;
    vector<vector<double>> mat2;
    vector<double> rowVector2(v);

    cout << "This is the second matrix: " << endl;
    while (myfile.good()) {
        mat2.push_back(rowVector2);
        for (int col = 0; col < v; col++) {
            myfile >> mat2[row2][col];
        }
        row2++;
        if (row2 >= u) {
            break;
        }
    }
    cout << mat2 << '\n' << endl;

    vector<vector<double>> product;
    product = dotprod(mat1, mat2, x, y, v);

    cout << product << '\n' << endl;
*/

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

/*
    // Encrypting matrix of integers.
    cout << Encrypt(m, p, r, L, c, w, x, v, product) << endl;
*/
