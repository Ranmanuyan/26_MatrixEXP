#include <iterator>
#include <string>
#include <ctime>
#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <chrono>


double cdfSumExp(std::vector<double> lambdaVec, double x){
  int n = lambdaVec.size();
  
  double prodLam = 1;
  for (int i = 0; i < n; ++i){
    prodLam *=  lambdaVec[i];    
  }


  double jNum = 0;
  for (int j = 0; j < n; ++j){
    double expNum = (exp(-1 * x * lambdaVec[j]))/(lambdaVec[j]);

    double lambDifProd = 1;
    for (int k = 0; k < n; ++k){
      if(k != j){
        lambDifProd *= (lambdaVec[k] - lambdaVec[j]);
      }
    }
    jNum += expNum/lambDifProd;
  }

  double cdfVal = 1- (jNum * prodLam);

  return cdfVal; 
}


int main(int argc, char const *argv[]){
  double x = 0;
  std::vector<double> lambdaVec {
    0.4, 0.7,0.3,0.7,0.3
  };
  
  double c = cdfSumExp(lambdaVec, x);
  
  std::cout << c <<std::endl;
  
  return 0;
}