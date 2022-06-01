
/* 
 * File:   expv.cpp
 * Author: Lei Liu
 * based on matlab expv, modify the krylov code to improve the accurancy
 *
 */


// NOTE ::
// * This read_data function, the last row is useless, i.e., in the data file, you need put an extra row at the end
// * then the collected data are correct, so in the data file (n+1) * n
//  not because of phase-type, just because of the read data readMatrix function


#include <Eigen/unsupported/Eigen/MatrixFunctions>
#include <Eigen/Core>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>       /* exp */
using namespace std;
using namespace Eigen;



// input: t, A, v, tol, m
// tol = 1e-7, m = 30, v = [1,1,1,1,1,1]^T

void expv(MatrixXd A,  double t, VectorXd &w){
    const double PI = std::atan(1.0)*4;
    double tol = 1e-4;
    double anorm = A.lpNorm<Eigen::Infinity>();
    int m = 30;
    int mb = m;
    int mxrej = 10;
    int n = A.cols();
    VectorXd v = VectorXd::Ones(n);

    int k1 = 2; 
    double xm = 1/double(m); 
    double normv = v.norm(); 
    double beta = normv;

    double fact = pow( (double((m+1)/exp(1.0))) , (m+1)) * sqrt(2*PI*(m+1));
    double t_new = (1/anorm)*pow(((fact*tol)/(4*beta*anorm)),xm);
    double s = pow(10,(floor(log10(t_new))-1)); 
    t_new = ceil(t_new/s)*s; 
    // int sgn = 1; 
    int nstep = 0;

    w = v;
    double hump = normv;
    double t_now =0;
    double avnorm = 0;
    double err_loc = tol;
    double phi1 = 0;
    double phi2 = 0;
    double delta = 1.2; 
    double gamma = 0.9;
    int mx =0;
    double s_error = 0;
    
    while (t_now < t){
        MatrixXd F;
        nstep = nstep + 1;
        double t_step = min( t-t_now,t_new);
        MatrixXd V = MatrixXd::Zero(n,m+1); 
        MatrixXd H = MatrixXd::Zero(m+2,m+2);
        V.col(0) = (1/beta)*w; // the first column
        for (int j = 0; j < m; ++j){
            VectorXd p = A * V.col(j);
            for (int i = 0; i <= j; ++i){  // check here, should be <= or < ?
                H(i,j) = V.col(i).dot(p);
                p = p - H(i,j)*V.col(i);
            }
            double s = p.norm();
            if (s < tol){
                k1 = 0;
                mb = j;
                t_step = t-t_now;
                break;
            }
            H(j+1,j) = s;
            V.col(j+1) = (1/s)*p;
        }
        if (k1 == 0){
            H(m+1,m) = 1;
            avnorm = (A*(V.col(m))).norm();
        }
        int ireject = 0;
        while (ireject <= mxrej){
            mx = mb + k1;
            F = (t_step * H.block(0,0,mx,mx)).exp();
            if (k1 == 0){
                err_loc = tol; 
                break;
            }
            else{
                phi1 = abs( beta*F(m,0) );
                phi2 = abs( beta*F(m+1,0) * avnorm );
            }
            if (phi1 > 10*phi2){
               err_loc = phi2;
               xm = 1/double(m);
            }
            else if (phi1 > phi2){
               err_loc = (phi1*phi2)/(phi1-phi2);
               xm = 1/double(m); 
           }else{
                err_loc = phi1;
                xm = 1/double(m-1);
           }
           if (err_loc <= delta * t_step*tol){
                break;
            }else{
                t_step = gamma * t_step * pow(t_step*tol/err_loc,xm);
                s = pow(10,(floor(log10(t_step))-1));
                t_step = ceil(t_step/s) * s;
                ireject = ireject + 1;
            }
    }
    mx = mb + max( 0,k1-1 );
    w = V.block(0,0,V.rows(),mx)*(beta*F.block(0,0,mx,1));
    beta = w.norm();
    hump = max(hump,beta);

    t_now = t_now + t_step;
    t_new = gamma * t_step * pow((t_step*tol/err_loc),xm);
    s = pow(10,(floor(log10(t_new))-1)); 
    t_new = ceil(t_new/s) * s;

    // err_loc = max(err_loc,rndoff);
    s_error = s_error + err_loc;
    } // end while 

double err = s_error;
hump = hump / normv;
}




















// ====================================================================================
/**
 * Arnoldi method to compute Krylov approximation of a matrix
 * @param A nxn matrix to approximate
 * @param v nx1 initial vector
 * @param k number of Krylov steps (size of resulting basis)
 * @param V output matrix (n x k) of orthogonal vectors
 * @param H output matrix (k+1 x k) containing Krylov approximation of A
 */
// template<class DerivedMatrix, class DerivedVector>

void arnoldi(MatrixXd A,  VectorXd v, int k, MatrixXd &V, MatrixXd &H){
    int n = A.cols();
    VectorXd vt(n);
    V.col(0) = v/v.norm();
    for (unsigned int m=0; m<k; m++) {
        vt.noalias() = A*V.col(m); // noalias can speed up the matrix product?
        for (unsigned int j=0; j<m+1; j++) {
            H(j,m) = vt.dot(V.col(j));
            vt = vt - H(j,m)*V.col(j);
        }
        H(m+1,m) = vt.norm();
        if (m != k-1)
            V.col(m+1) = vt/H(m+1,m);
    }
}



double cdfCal(MatrixXd Qstar, VectorXd a, double x1){
    MatrixXd mat = (x1 * Qstar).exp();
    double cdf = 1-((a.transpose() * mat).sum());
    mat.resize(0,0);
    return cdf;
}

#define MAXBUFSIZE  ((int) 1e6)

MatrixXd readMatrix(const char *filename){
    int cols = 0, rows = 0;
    double buff[MAXBUFSIZE];

    // Read numbers from file into buffer.
    ifstream infile;
    infile.open(filename);
    while (! infile.eof())
        {
        string line;
        getline(infile, line);

        int temp_cols = 0;
        stringstream stream(line);
        while(! stream.eof())
            stream >> buff[cols*rows+temp_cols++];

        if (temp_cols == 0)
            continue;

        if (cols == 0)
            cols = temp_cols;

        rows++;
        }

    infile.close();

    rows--;

    // Populate matrix with numbers.
    MatrixXd result(rows,cols);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < cols; j++)
            result(i,j) = buff[ cols*i+j ];

    return result;
};


int main(int argc, char const *argv[]){

    // read matrix data from file
    MatrixXd A = readMatrix("Mat31.txt");
    int n = A.cols();

    VectorXd a(n);
    a.setZero();
    a[0]= 1.0;
    for (int i = 1; i < n; ++i) {
        a[i] = 0;
    }

    int k = n/2;
    VectorXd v(n);
    for (int i = 0; i < n; ++i) {
        v[i] = 1;
    }

    VectorXd e1(k);
    e1.setZero();
    e1[0]= 1.0;
    for (int i = 1; i < k; ++i) {
        e1[i] = 0;
    }



   // The EIGEN Version to calculate the cdf of 100.
   std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();
   double cdf = cdfCal(A, a, 353);
   std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
   std::cout << "Time Costs EXPM = " << std::chrono::duration_cast<std::chrono::microseconds>(end1 - begin1).count() << "[microseconds]" << std::endl;
   std::cout << cdf << std::endl;

    // The expv version
   VectorXd w1 = VectorXd::Ones(n);
   // VectorXd &w = w1;
   std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
   expv( A,  353, w1);
   std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
   std::cout << 1 - w1(0) << std::endl;
   std::cout << "Time Costs Krylov= " << std::chrono::duration_cast<std::chrono::microseconds>(end2 - begin2).count() << "[microseconds]" << std::endl;





   //  /*
   //  *
   //  *   The Krylov subspace version, with sub_matrix size(30*30), 
   //  *   it's about half of the original matrix siez 57
   //  *
   //  *
   //  */
   

   // // * @param V output matrix (n x k) of orthogonal vectors
   // //  * @param H output matrix (k+1 x k) containing Krylov approximation of A

   // std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
   // //Please control the non-used elements: should be 0, H is a upper Hessenberg
   // MatrixXd H = MatrixXd::Zero(k+1, k);
   // MatrixXd V = MatrixXd::Zero(n, k);
   // arnoldi(A, v, k,  V,  H);
   // // here we have (k+1 * k), H_(k+1, k)==?
   
   // MatrixXd H1 = H.topRows(k);
   // double cdf2 =  1- (v.norm() * V *(((420*H1).exp())* e1))(0);
   // std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
   // std::cout << "Time Costs Krylov= " << std::chrono::duration_cast<std::chrono::microseconds>(end2 - begin2).count() << "[microseconds]" << std::endl;
   // std::cout << cdf2 << std::endl;

    return 0;
}


// EIGEN optimizer from zhihu:   -mavx -mfma
// g++ -O0 -std=c++11 -g  expv.cpp -o expv

// Compile with this one!!
// g++ -O3 -std=c++11 -g  -march=native expv.cpp -o expv