

#include <Eigen/unsupported/Eigen/MatrixFunctions>
#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#include <algorithm> 
#include <math.h>       /* pow */
#include <vector>

using namespace Eigen;

// https://stackoverflow.com/questions/43151853/eigen-balancing-matrix-for-eigenvalue


void balance_matrix(Eigen::MatrixXd &A, Eigen::MatrixXd &Aprime, Eigen::MatrixXd &D) {
    // https://arxiv.org/pdf/1401.5766.pdf (Algorithm #3)
    const int p = 2;
    double beta = 2; // Radix base (2?)
    Aprime = A;
    D = Eigen::MatrixXd::Identity(A.rows(), A.cols());
    bool converged = false;
    do {
        converged = true;
        for (Eigen::Index i = 0; i < A.rows(); ++i) {
            double c = Aprime.col(i).lpNorm<p>();
            double r = Aprime.row(i).lpNorm<p>();
            double s = pow(c, p) + pow(r, p);
            double f = 1;
            while (c < r / beta) {
                c *= beta;
                r /= beta;
                f *= beta;
            }
            while (c >= r*beta) {
                c /= beta;
                r *= beta;
                f /= beta;
            }
            if (pow(c, p) + pow(r, p) < 0.95*s) {
                converged = false;
                D(i, i) *= f;
                Aprime.col(i) *= f;
                Aprime.row(i) /= f;
            }
        }
    } while (!converged);
}

int main(int argc, char const *argv[]){

    MatrixXd A(29,29);
    A << -0.48,0.48,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
         0,-0.68,0.41,0.27,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
         0,0,-0.27,0,0.27,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.41,0.41,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.8,0.28,0.52,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.52,0,0.52,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.28,0.28,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.7,0.53,0.17,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.17,0,0.17,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.53,0.53,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.96,0.87,0.09,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.09,0,0.09,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.87,0.87,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.68,0.61,0.07,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.07,0,0.07,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.61,0.61,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.34,0.85,0.49,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.49,0,0.49,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.85,0.85,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.95,0.92,0.03,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.03,0,0.03,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.92,0.92,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.99,0.8,0.19,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.19,0,0.19,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.8,0.8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1.08,0.9,0.18,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.18,0,0.18,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.9,0.9,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-0.16;
    MatrixXd Aprime = A;
    MatrixXd D = Eigen::MatrixXd::Identity(A.rows(), A.cols());

    // std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();

    balance_matrix(A, Aprime, D);
    // std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
    // std::cout << "Time Costs = " << std::chrono::duration_cast<std::chrono::microseconds>(end1 - begin1).count() << "[microseconds]" << std::endl;

    std::cout<< D << std::endl;

    std::cout<< "=============" << std::endl;

    std::cout<< Aprime << std::endl;

    std::cout<< "=============" << std::endl;

    std::cout<< D.inverse()*A*D << std::endl;



    return 0;
}

