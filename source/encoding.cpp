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

vector<vector<double>> matrix_transpose(vector<vector<double>> mat1) {

    vector<vector<double>> product_trans(mat1[0].size(), vector<double>(mat1.size()));

    for (size_t i = 0; i < mat1.size(); ++i)
        for (size_t j = 0; j < mat1[0].size(); ++j)
            product_trans[j][i] = mat1[i][j];

    return product_trans;

}

vector<vector<double>> dotprod(vector<vector<double>> mat1, vector<vector<double>> mat2, int x, int y) {

    vector<vector<double>> mult;

    // Dot Product.
    mult.resize(y);
    for (int i = 0; i < mult.size(); i++) {
        mult[i].resize(y);
    }

//    cout << "Dot Product of Matrices: " << endl;
    for (int i = 0; i < y; i++) {
        for (int j = 0; j < y; j++) {
            mult[i][j] = 0;
            for (int k = 0; k < x; k++) {
                mult[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }

    return mult;

}

// Inverse function.
vector<vector<double>> Inv(int y, vector<vector<double>> mult) {

    double t, div;

    mult.resize(y);
    for (int i = 0; i < mult.size(); i++) {
        mult[i].resize(2 * y);
    }

    for(int i = 0; i < y; i++) {
        for(int j = y; j < 2 * y; j++) {
            if (i == j - y) {
                mult[i][j] = 1;
            }
            else {
                mult[i][j] = 0;
            }
        }
    }

    for(int i = 0; i < y; i++)
    {
        t = mult[i][i];
        for(int j = i; j < 2 * y; j++)
            mult[i][j]=mult[i][j]/t;
        for(int j = 0; j < y; j++)
        {
            if(i!=j)
            {
                t=mult[j][i];
                for(int k = 0; k < 2 * y; k++)
                    mult[j][k] = mult[j][k] - (t * mult[i][k]);
            }
            div = mult[0][0];
        }
    }

    cout<<"Inverse matrix of Xtrans_X is: ";
    for(int i = 0; i < y; i++) {
        for(int j = y; j < 2 * y; j++) {
            mult[i][j] = mult[i][j] / div;
        }
        for (int j = 0; j < 2 * y; j++) {
            // Bring forward the inverse matrix.
            mult[i][j] = mult[i][j+y];
            if (j+y >= 2 * y) {
                break;
            }
        }
    }

    // Resize the matrix to keep only the inverse matrix.
    mult.resize(y);
    for (int i = 0; i < mult.size(); i++) {
        mult[i].resize(y);
    }

    return mult;

}

vector<vector<ZZ>> Encrypt(long m, long p, long r, long L, long c, long w, int rows, int cols, vector<vector<double>> mat, vector<vector<double>> dec, int phim) {
//ZZX Encrypt(long m, long p, long r, long L, long c, long w, int rows, int cols, vector<vector<double>> mat, vector<vector<double>> dec, int phim) {

    FHEcontext context(m, p, r);
    buildModChain(context, L, c);
    FHESecKey secretKey(context);
    FHEPubKey &publicKey = secretKey;
    secretKey.GenSecKey(w);

    vector<vector<ZZX>> matrix;

    auto begin_encrypt = Clock::now();

    // Conversion of fractions (int part + dec part) to ZZX.
    matrix = frac_to_ZZX(rows, cols, mat, dec, phim);

    cout << matrix << '\n' << endl;

    // Encrypt for all ZZX in mat
    Ctxt enc(publicKey);
    vector<vector<Ctxt>> ctxt_mat;

    for (int i = 0; i < rows; i++) {
        vector<Ctxt> temp_ctxt;
        for (int j = 0; j < cols; j++) {
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

    //vector<vector<double>> mat_ans;
    vector<vector<ZZ>> mat_ans;
    //double temp_store;
    ZZX temp_store_zzx;


    ZZX pptxt;
    ZZ temp_storee;
    vector<ZZ> temp_store;
    vector<vector<ZZ>> ans;
    vector<ZZ> temp1_store;
    ZZ temp1_storee;
    ZZ final;

    int counter = 0;

    for (int i = 0; i < rows; i++) {
        //vector<ZZ> mat_temp;
        //vector<double> mat_temp;
        //vector<Ctxt> temp_ctxt;
        for (int j = 0; j < cols; j++) {
            secretKey.Decrypt(temp_store_zzx, ctxt_mat[i][j]);
            //decode temp_store_zzx into temp_store
            //decode(temp_store_zzx, rows, cols, mat, phim);
            cout << "This is zzx before loop: " << temp_store_zzx << endl;
//            for (int k = 0; k < temp_store_zzx.size(); k++) {

            for (int k = 0; k < phim; k++) {
                counter++;

                // Fractional part.
//                if (k > temp_store_zzx.size() / 2) {
                if (k > phim / 2) {
                    // for every element
                    if (temp_store_zzx[k] == coeff(temp_store_zzx, k)) {
                        // -x^(n-i) --> x^i
                        // SetCoeff(pptxt, -1 * (i - phim), -1 * coeff(pptxt, i)); --> **pptxt will be empty.
                        SetCoeff(temp_store_zzx, -(k - phim), -coeff(temp_store_zzx, k));
                    }
                    // Exponent is negative i. x^(-i) --> 1 / x^i.
                    double power = 1 / pow(2, k);
                    // temp1_storee = pptxt[i] * power;
                    temp1_storee += temp_store_zzx[k] * power;
//                    if (temp_store_zzx[i] == []) {
//                        break;
//                    }
                }
                    // Integer part.
                else {
                    temp_storee += temp_store_zzx[k] * pow(2, i);
//                    if (temp_store_zzx[i] == 0) {
//                        break;
//                    }
                }
                add(final, temp1_storee, temp_storee);
                //temp_store.push_back(final);
            }
            cout << "This is counter: " << counter << endl;
            //ans.push_back(temp_store);
            temp_store.push_back(final);
        }
        ans.push_back(temp_store);
    }












//    return temp_store_zzx;

//    cout << decode << endl;

/*
            mat_temp.push_back(finally);
        }
        mat_ans.push_back(mat_temp);
    }
*/
    cout << "-------------------- Decryption --------------------" << endl;
    cout << "Plaintext:  " << endl;

//    return mat_ans;
    return ans;



}

/*
//double decode(ZZX temp_store_zzx, int rows, int cols, vector<vector<double>> mat, int phim) {
//ZZX decode(ZZX temp_store_zzx, int rows, int cols, vector<vector<double>> mat, int phim) {
ZZ decode(ZZX temp_store_zzx, int rows, int cols, vector<vector<double>> mat, int phim) {

    ZZX pptxt;
    ZZ temp_storee;
    vector<ZZ> temp_store;
    vector<vector<ZZ>> what;

//    cout << "After decoding: " << endl;
//    for (int i = 0; i < phim/2; i++) {
//        temp_storee = temp_store_zzx[i] *= pow(2, i);
//        temp_store.push_back(temp_storee);
//        if (temp_store_zzx[i] == 0) {
//            break;
//        }
//    what.push_back(temp_store);
//    break;
//    }
//    cout << "print" << endl;
//    // Code cant exit the above loop. !!
//    for (int i = 0; i >= phim/2; i++) {
//        if (temp_store_zzx[i] == coeff(temp_store_zzx, i)) {
//            // change the sign of the term and minus power n.
//            SetCoeff(pptxt, i - phim, -1 * coeff(temp_store_zzx, i));
//        }
//    }

    vector<ZZ> temp1_store;
    vector<vector<ZZ>> whut;
    ZZ temp1_storee;
    ZZ finally;

    cout << "After decoding: " << endl;
    for (int i = 0; i < phim; i++) {
        // Fractional part.
        if (i > phim/2) {
            // DOESNT RUN AFTER HERE. !!
            // for every element
            if (temp_store_zzx[i] == coeff(temp_store_zzx, i)) {
                // change the sign of the term and minus power n.
                SetCoeff(pptxt, i - phim, -1 * coeff(temp_store_zzx, i));
            }
            temp1_storee = pptxt[i] * pow(2, i);
            //temp1_store.push_back(temp1_storee);
        }
        //whut.push_back(temp1_store);
        // Integer part.
        else {
            temp_storee = temp_store_zzx[i] * pow(2, i);
            finally = temp1_storee + temp_storee;
            temp_store.push_back(finally);
            if (temp_store_zzx[i] == 0) {
                break;
            }
            what.push_back(temp_store);
            break;
        }
    }

//    return temp_store_zzx;
    return finally;


}
*/


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
        }
    }

    return ptxt;

}

// Adding the int part tgt with the frac part.
//vector<vector<ZZX>> frac_to_ZZX(int rows, int cols, vector<vector<ZZX>> mat_int, vector<vector<ZZX>> binary, ZZX msg1, ZZX msg2) {
vector<vector<ZZX>> frac_to_ZZX(int rows, int cols, vector<vector<double>> mat, vector<vector<double>> dec, int phim) {

    ZZX msg1;
    ZZX msg2;
    ZZX int_and_frac;
    vector<vector<ZZX>> polyn;

    cout << "fraction to ZZX: " << endl;
    for (int i = 0; i < rows; i++) {
        vector<ZZX> tempp;
        for (int j = 0; j < cols; j++) {
            int z1 = mat[i][j];
            // Int part --> ZZX.
            msg1 = encode(z1);
            long double z2 = dec[i][j];
            // Fractional part --> ZZX.
            msg2 = frac_encoder(z2, cols, phim);
            // Adding the int and fractional parts tgt.
            int_and_frac = msg1 + msg2;
            tempp.push_back(int_and_frac);
        }
        polyn.push_back(tempp);
    }

    return polyn;

}












/*
vector<vector<double>> Y() {

    int x, y;

    ifstream myfile;
    myfile.open("/home/karis/CLionProjects/HElib-basic/matrix.txt");

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

    cout << "This is the matrix Y: ";
    for (int i = 0; i < x; i++) {
        mat1[i][0] = 1;
        for (int j = 1; j < y; j++) {
            myfile >> mat1[i][j];
        }
        myfile >> mat2[i][0];
    }

    return mat2;

}

vector<vector<double>> Xtrans_X() {

    int x, y;

    ifstream myfile;
    myfile.open("/home/karis/CLionProjects/HElib-basic/matrix.txt");

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

    cout << "This is the matrix X: " << mat1 << endl;
//    cout << "This is the matrix Y: " << mat2 << endl;

    vector<vector<double>> mat1_trans;
    mat1_trans = matrix_transpose(mat1);
    cout << "This is X transpose: " << mat1_trans << endl;

    vector<vector<double>> product;
    product = dotprod(mat1_trans, mat1, x, y);
    cout << "This is Xtrans_X: ";
    return product;

}

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