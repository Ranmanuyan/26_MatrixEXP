
/* 
 * File:   newArnoldi.cpp
 * Author: Lei Liu
 * the arnoldi code referred from github of Robert Gantner Created on November 25, 2011
 *
 */


// NOTE ::
// * This read_data function, the last row is useless, i.e., in the data file, you need put an extra row at the end
// * then the collected data are correct
//  not phase-type, just because of the read data readMatrix function


#include <Eigen/unsupported/Eigen/MatrixFunctions>
#include <Eigen/Core>
#include <vector>
#include <iostream>
#include <fstream>
using namespace std;
using namespace Eigen;


/**
 * Arnoldi method to compute Krylov approximation of a matrix
 * @param A nxn matrix to approximate
 * @param v nx1 initial vector
 * @param k number of Krylov steps (size of resulting basis)
 * @param V output matrix (n x k) of orthogonal vectors
 * @param H output matrix (k+1 x k) containing Krylov approximation of A
 */
// template<class DerivedMatrix, class DerivedVector>

void arnoldi( MatrixXd A,  VectorXd v, int k, MatrixXd &V, MatrixXd &H){
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

// for test to see the matrix
// MatrixXd matCal(MatrixXd Qstar, VectorXd a, double x1){
//     MatrixXd mat = (x1 * Qstar).exp();
//     return mat;
// }






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
    MatrixXd A = readMatrix("Mat169.txt");
    int n = A.cols();
    int k = n/2;

    VectorXd a(n);
    a.setZero();
    a[0]= 1.0;
    for (int i = 1; i < n; ++i) {
        a[i] = 0;
    }
    
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
   double cdf = cdfCal(A, a, 450);
   std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
   std::cout << "Time Costs EXPM = " << std::chrono::duration_cast<std::chrono::microseconds>(end1 - begin1).count() << "[microseconds]" << std::endl;
   std::cout << cdf << std::endl;


    /*
    *
    *   The Krylov subspace version, with sub_matrix size(30*30), 
    *   it's about half of the original matrix siez 57
    *
    *
    */
   

   // * @param V output matrix (n x k) of orthogonal vectors
   //  * @param H output matrix (k+1 x k) containing Krylov approximation of A

   std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
   //Please control the non-used elements: should be 0, H is a upper Hessenberg
   MatrixXd H = MatrixXd::Zero(k+1, k);
   MatrixXd V = MatrixXd::Zero(n, k);
   arnoldi(A, v, k,  V,  H);
   // here we have (k+1 * k), H_(k+1, k)==?
   
   MatrixXd H1 = H.topRows(k);

   double cdf2 =  1- (v.norm() * V *(((450*H1).exp())* e1))(0);
   std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
   std::cout << "Time Costs Krylov= " << std::chrono::duration_cast<std::chrono::microseconds>(end2 - begin2).count() << "[microseconds]" << std::endl;
   std::cout << cdf2 << std::endl;

    //  first mat*t
   // std::chrono::steady_clock::time_point begin3 = std::chrono::steady_clock::now();
   // MatrixXd Hm(31, 30);
   // MatrixXd Vm(57, 30);
   // arnoldi(100*A, v, 30,  Vm,  Hm);
   // MatrixXd H1m = Hm.topRows(30);
   // double cdf3 =  1- (v.norm() * Vm *(((H1m).exp())* e1))(0);
   // std::chrono::steady_clock::time_point end3 = std::chrono::steady_clock::now();
   // std::cout << "Time Costs = " << std::chrono::duration_cast<std::chrono::microseconds>(end3 - begin3).count() << "[microseconds]" << std::endl;
   // std::cout << cdf3 << std::endl;





    return 0;
}






   // double cdf2 =  1- (a.transpose() * v.norm() * V * ((100*H1).exp())*e1);
   // std::cout << cdf2 << std::endl;
   // std::cout << V  << std::endl;
   // std::cout << "------------------"  << std::endl;
   // std::cout << (v.norm())<< std::endl;
   // std::cout << "------------------"  << std::endl;    
    // // MatrixXd a2 = matCal(A, a, 100);

    


    // MatrixXd Vm = V.topRows(10);

    // MatrixXd T1 = (V.transpose()*V).inverse();
    // MatrixXd T2 = (Vm*Vm.transpose());
    // MatrixXd T3 = Vm*((Vm.transpose()*Vm).inverse()) * (Vm.transpose());
    // // std::cout << (V*((V.transpose()*V).inverse()) * (V.transpose()) * a2 * v) << std::endl;
    // std::cout << Vm << std::endl;

    // std::cout << (a2 * v) << std::endl;





    // std::cout << V * ((V.transpose()*V).inverse()) * V.transpose() << std::endl;
    // std::cout<< "-----------------------------" <<std::endl;

    // std::cout<< v.norm()* V*(V.transpose()*V).inverse()* V.transpose() * a2 *V * e1 <<std::endl;
    // std::cout<< "-----------------------------" <<std::endl;
    // std::cout<< V.transpose() <<std::endl;
    // std::cout<< "-----------------------------" <<std::endl;
    // std::cout<< V.transpose() * V <<std::endl;




    // std::cout<< H1.exp() <<std::endl;
    // std::cout<< "-----------------------------" <<std::endl;
    // std::cout << v.norm() << std::endl;
    // std::cout << (100*H1).exp() << std::endl;
    // std::cout << "-----------------------------" << std::endl;
    // MatrixXd x =  v.norm() * V * ((100*H1).exp());
    // std::cout << x << std::endl;
    // std::cout << "-----------------------------" << std::endl;
    // MatrixXd x1 =  V.transpose() * (A.exp()) * V;
    // std::cout << "-----------------------------" << std::endl;
    // std::cout << x1 << std::endl;
    // std::cout << "-----------------------------" << std::endl;
    // std::cout << v.norm() * V << std::endl;


    // double cdf2 = 1 - cdfCal(A, a, 100);
    // MatrixXd a2 = matCal(A, a, 100);

    // std::cout << a2 << std::endl;





// EIGEN optimizer from zhihu:   -mavx -mfma
// g++ -O0 -std=c++11 -g  newArondi.cpp -o newArnoldi
// g++ -O3 -std=c++11 -g  -march=native newArondi.cpp -o newArnoldi