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

#include <string>
#include <sstream>
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

    cout << X() << endl;
    cout << Y() << endl;
/*    int x, y, u, v;

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

/*
    // Encrypting matrix of integers.
    cout << Encrypt(m, p, r, L, c, w, x, v, product) << endl;
*/
    // Fractional Encoder.

    cout << "-------------------- Encryption --------------------" << endl;
    auto begin_encrypt = Clock::now();

    int rows, cols;

    ifstream infile;
    infile.open("/home/karis/CLionProjects/HElib-basic/dec.txt");

    // Tests if the file opens successfully.
    if (!infile.is_open()) {
        cout << "File failed to open!" << endl;
    }

    infile >> rows;
    cout << rows << endl; // no. of rows
    infile >> cols;
    cout << cols << endl; // no. of cols

    int row_count = 0;
    vector<vector<double>> mat; // 2d array as a vector of vectors
    vector<double> rowvec(cols);

    // Storing the decimal matrix in mat.
    cout << "This is the decimal matrix: " << endl;
    while (infile.good()) {
        mat.push_back(rowvec);
        for (int col = 0; col < cols; col++) {
            infile >> mat[row_count][col];
        }
        row_count++;
        if (row_count >= rows) {
            break;
        }
    }
    cout << mat << '\n' << endl;

    vector<vector<double>> dec;

    dec = mat;

    // The remaining decimal value: initial value - integer.
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            dec[i][j] = dec[i][j] - floor(dec[i][j]);
        }
    }

    cout << "This is the remaining dec value: " << '\n' << dec << '\n' << endl;

    // Rounding down the elements in mat to integers.
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            mat[i][j] = floor(mat[i][j]);
            cout << mat[i][j] << " ";
        }
        cout << '\n' << endl;
    }

    // Matrix rounded down to integers.
    cout << "Elements in matrix rounded down to int: " << '\n' << mat << '\n' << endl;

    // Conversion of integral part into binary/polynomial.
    cout << int_to_ZZX(rows, cols, mat) << '\n' << endl;

    // To get phim.
    int phim = to_int(context.zMStar.getPhiM());
    cout << phim << '\n' << endl;

    double z = 0.875;
    cout << frac_encoder(z, cols, phim)  << '\n' << endl;

    // Conversion of fractional part into binary.
    cout << frac_to_binary(rows, cols, dec, phim) << '\n' << endl;

    // Conversion of frac part to ZZX.
    // Adding the int part and fractional part together.
    //vector<vector<ZZX>> a = int_to_ZZX(rows, cols, mat);
    //vector<vector<ZZX>> b = frac_to_binary(rows, cols, dec, phim);
    //cout << frac_to_ZZX(m, p, r, a, b) << endl;















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


}