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

vector<vector<double>> dotprod(vector<vector<double>> mat1, vector<vector<double>> mat2, int x, int y, int z) {

    int mat1row = mat1.size();
    int mat1col = mat1[0].size();

    vector<vector<double>> mult;
    mult.resize(x, vector<double>(y));

//    cout << "Dot Product of Matrices: " << endl;
    for (int i = 0; i < x; i++) {
        for (int j = 0; j < y; j++) {
            mult[i][j] = 0;
            for (int k = 0; k < z; k++) {
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

vector<vector<Ctxt>> mat_mat_mult() {



}

vector<vector<Ctxt>> mat_vec_mult() {




}