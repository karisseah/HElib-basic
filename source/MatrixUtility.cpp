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

vector<vector<Ctxt>> mat_mat_mult(vector<vector<Ctxt>> mat1, vector<vector<Ctxt>> mat2) {

    int mat1row = mat1.size();

    int mat2row = mat2.size();
    int mat2col = mat2[0].size();

    // Once resized, errors.
    vector<vector<Ctxt>> mat;
    mat.resize(mat1row);
    for (int i = 0; i < mat1row; i++) {
        mat[i].resize(mat2col);
    }

    cout << "Ciphertext after multiplication:" << endl;
    for (int i = 0; i < mat1row; i++) {
        cout << "1st loop i = " << i << endl;
        for (int j = 0; j < mat2col; j++) {
            cout << "1st loop j = " << j << endl;
            for (int k = 0; k < mat2row; k++) {
                cout << "1st loop k = " << k << endl;
//                mat1[i][k] *= (mat2[k][j]);
                mat1[i][k].multiplyBy(mat2[k][j]);
                cout << mat1[i][k] << endl;
                mat[i][j] = mat1[i][k];
//                mat1[i][k].reLinearize(-1);
                cout << "finish mult1" << endl;
            }
        }
    }

//    for (int i = 0; i < mat1row; i++) {
//        for (int j = 0; j < mat2col; j++) {
//            mat[i][j] = mat1[i][j];
//        }
//    }

    return mat2;

}

vector<vector<Ctxt>> mat_vec_mult(vector<vector<Ctxt>> mat1, vector<vector<Ctxt>> mat2) {

    int mat1row = mat1.size();

    int mat2row = mat2.size();
    int mat2col = 1;

    for (int i = 0; i < mat1row; i++) {
        cout << "2nd loop i = " << i << endl;
        for (int j = 0; j < mat2col; j++) {
            cout << "2nd loop j = " << j << endl;
            for (int k = 0; k < mat2row; k++) {
                cout << "2nd loop k = " << k << endl;
//                mat1[i][k] *= (mat2[k][j]);
                mat1[i][k].multiplyBy(mat2[k][j]);
//                mat1[i][k].reLinearize(-1);
                cout << "finish mult2" << endl;
            }
        }
    }

    return mat1;

}