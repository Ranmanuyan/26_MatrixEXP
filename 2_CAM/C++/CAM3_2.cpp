// in this fragment, we tested  Code Fragment 3.2
// 21th Feb 2022

#include <fstream>
#include <Eigen/unsupported/Eigen/MatrixFunctions>
#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#include <algorithm> 
#include <math.h>       /* pow */
#include <vector>

using namespace std;
using namespace Eigen;


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


void calAlpha(MatrixXd A, std::vector<int> mList, std::vector<double> thetaList, int &sNum, int &numM){
  for (int p = 2; p < 9; ++p){
    // cal the alpha(A) first
    double d1 = pow(A.pow(p).lpNorm<1>(), (1.0/p));
    double d2 = pow(A.pow(p+1).lpNorm<1>(), (1.0/(p+1)));
    double alphaA = std::max(d1, d2);
    int index =0;

    for (int i = 0; i < mList.size(); ++i){
      if(mList[i] * ceil(alphaA/thetaList[i]) <  numM * ceil(alphaA/thetaList[index])){
        index = i;
        numM = mList[i];
        sNum = std::max(int(ceil(alphaA/thetaList[index])), 1);
       }
     }
    }
}

// function to calculate number s and number m
void sANDm(int &sNum, int &numM, MatrixXd A){
  // These are prepared based on the precision, no need to change.
  std::vector<int> mList {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55};
  std::vector<double> thetaList  {0.13, 1.0, 2.2, 3.6, 4.9, 6.3, 7.7, 9.1, 11.0, 12.0, 13.0};
  
  double normA = A.lpNorm<1>();
  int index = 0;

  // this number 41.6 is predefined
  if (normA <= 41.6){
    for (int i = 1; i < mList.size(); ++i){
      if(mList[i] * int(ceil(normA/thetaList[i])) <  numM * int(ceil(normA/thetaList[index]))){
        index = i;
        numM = mList[i];
      }
    }
    sNum = int(ceil(normA/thetaList[index]));
  }
  else{ // NOT validate yet
    calAlpha( A,  mList, thetaList,  sNum,  numM);
  }
}

void CAM_algo(MatrixXd &A, VectorXd &v, double t, VectorXd &F){
   // balance and get Aprime, D, D^-1,
   MatrixXd Aprime = A;
   MatrixXd D = Eigen::MatrixXd::Identity(A.rows(), A.cols());
   balance_matrix(A, Aprime, D);
   
   // compare the norm, and get B
   if(Aprime.lpNorm<1>()<A.lpNorm<1>()){
    A = Aprime;
    v = D.inverse()*v;
   }
   double traceVal = A.trace()/A.rows();
   A = A - traceVal*Eigen::MatrixXd::Identity(A.rows(), A.cols());

   // calculate s, m
   int m= 5;
   int s = 1;
   sANDm(s, m, t*A);

   // begin F algorithm
      
   F = v;
   double ita = exp(t*traceVal/s);
   for (int i = 1; i < s+1; ++i){
     double c1 = v.lpNorm<Infinity>();
     for (int j = 1; j < m+1; ++j){
       v = (t*(A*v))/(s*j);
       double c2 = v.lpNorm<Infinity>();
       F = F+v;
       if(c1+c2 <= 6e-8 * F.lpNorm<Infinity>()){
        break;
       }
       c1 = c2;
     }
     F = ita * F;
     v = F;
   }
   F = D*F;

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

int main(){
  
  MatrixXd A = readMatrix("Mat31.txt");
  int n = A.cols();
  VectorXd v = VectorXd::Ones(n);
  double t = 420;
   

   std::chrono::steady_clock::time_point begin1 = std::chrono::steady_clock::now();
   MatrixXd final1 =  (A.exp().pow(t))*v;
   std::cout << 1- final1(0) <<std::endl;
   std::chrono::steady_clock::time_point end1 = std::chrono::steady_clock::now();
   std::cout << "Time Costs = " << std::chrono::duration_cast<std::chrono::microseconds>(end1 - begin1).count() << "[microseconds]" << std::endl;  

   std::chrono::steady_clock::time_point begin2 = std::chrono::steady_clock::now();
   VectorXd F = v;
   CAM_algo(A, v, t, F);
   std::cout << 1- F(0) <<std::endl;
   std::chrono::steady_clock::time_point end2 = std::chrono::steady_clock::now();
   std::cout << "Time Costs = " << std::chrono::duration_cast<std::chrono::microseconds>(end2 - begin2).count() << "[microseconds]" << std::endl;



}
// compile command
// g++ -O0 -std=c++11 -g  CAM3_2.cpp -o test1