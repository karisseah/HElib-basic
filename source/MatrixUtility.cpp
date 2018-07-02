//
// Created by karis on 6/29/18.
//

#include "MatrixUtility.h"

vector<vector<double>> matrix_transpose(vector<vector<double>> mat1) {

    vector<vector<double>> product_trans(mat1[0].size(), vector<double>(mat1.size()));

    for (size_t i = 0; i < mat1.size(); ++i)
        for (size_t j = 0; j < mat1[0].size(); ++j)
            product_trans[j][i] = mat1[i][j];

    return product_trans;

}

vector<vector<double>> dotprod(vector<vector<double>> mat1, vector<vector<double>> mat2) {

    int mat1row = mat1.size();
    int mat1col = mat1[0].size();

    int mat2row = mat2.size();
    int mat2col = mat2[0].size();

    vector<vector<double>> mult;
    mult.resize(mat1row, vector<double>(mat2col));

//    cout << "Dot Product of Matrices: " << endl;
    for (int i = 0; i < mat1row; i++) {
        for (int j = 0; j < mat2col; j++) {
            mult[i][j] = 0;
            for (int k = 0; k < mat2row; k++) {
                mult[i][j] += mat1[i][k] * mat2[k][j];
            }
        }
    }

    return mult;

}

// Inverse function.
vector<vector<double>> Inv(vector<vector<double>> mult) {

    int multsize = mult.size();

    double t, div;

    mult.resize(multsize);
    for (int i = 0; i < multsize; i++) {
        mult[i].resize(2 * multsize);
    }

    for(int i = 0; i < multsize; i++) {
        for(int j = multsize; j < 2 * multsize; j++) {
            if (i == j - multsize) {
                mult[i][j] = 1;
            }
            else {
                mult[i][j] = 0;
            }
        }
    }

    for(int i = 0; i < multsize; i++) {
        t = mult[i][i];
        for(int j = i; j < 2 * multsize; j++)
            mult[i][j] = mult[i][j] / t;
        for(int j = 0; j < multsize; j++) {
            if(i != j) {
                t = mult[j][i];
                for(int k = 0; k < 2 * multsize; k++)
                    mult[j][k] = mult[j][k] - (t * mult[i][k]);
            }
            div = mult[0][0];
        }
    }

    for(int i = 0; i < multsize; i++) {
        for(int j = multsize; j < 2 * multsize; j++) {
            mult[i][j] = mult[i][j] / div;
        }
        for (int j = 0; j < 2 * multsize; j++) {
            // Bring forward the inverse matrix.
            mult[i][j] = mult[i][j + multsize];
            if (j + multsize >= 2 * multsize) {
                break;
            }
        }
    }

    // Resize the matrix to keep only the inverse matrix.

    mult.resize(multsize);
    for (int i = 0; i < multsize; i++) {
        mult[i].resize(multsize);
    }

    return mult;

}

vector<vector<Ctxt>> mat_mat_mult(vector<vector<Ctxt>> mat1, vector<vector<Ctxt>> mat2, vector<vector<Ctxt>> mat3) {

    int mat1row = mat1.size();

    int mat2row = mat2.size();
    int mat2col = mat2[0].size();

    int mat3col = mat3[0].size();

    mat1[0][0].multiplyBy(mat2[0][0]);

//    cout << "Ciphertext after multiplication:" << endl;
//    for (int i = 0; i < mat1row; i++) {
//        cout << "1st loop i = " << i << endl;
//        for (int j = 0; j < mat2col; j++) {
//            cout << "1st loop j = " << j << endl;
//            for (int k = 0; k < mat2row; k++) {
//                cout << "1st loop k = " << k << endl;
//                mat1[i][k] *= (mat2[k][j]);
//                mat1[i][k].multiplyBy(mat2[k][j]);
//                mat1[i][k].reLinearize(-1);
//                //ctxt1_mat[i][j].multiplyBy2(ctxt2_mat[i][j], ctxt3_mat[i][j]);
//                cout << "finish mult1" << endl;
//            }
//        }
//    }

    for (int i = 0; i < mat1row; i++) {
        cout << "i = " << i << endl;
        for (int j = 0; j < mat3col; j++) {
            cout << "j = " << j << endl;
            for (int k = 0; k < mat2col; k++) {
                cout << "k = " << k << endl;
//                mat1[i][k].multiplyBy(mat3[k][j]);
                //ctxt1_mat[i][j].multiplyBy2(ctxt2_mat[i][j], ctxt3_mat[i][j]);
                cout << "finish mult2" << endl;
            }
        }
    }

    return mat1;

}