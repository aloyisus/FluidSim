#include <iostream>
#include <openvdb/math/ConjGradient.h>
#include "grid.h"
#include "maths.h"



bool CGSolver(const SymmBandMatrix& A, const SymmBandMatrix::VectorType& b, SymmBandMatrix::VectorType& p, double tol, int maxIter){
    
    p.fill(0);

    SymmBandMatrix::VectorType r = b;
    
    SymmBandMatrix::VectorType s(b.size(),0);
    
    SymmBandMatrix::VectorType z(b.size(),0);

    s = r;
    
    double sigma = r.dot(r);
    if (sigma == 0.0) {
        std::cout << "CG early exit: RHS is zero\n";
        return true;
    }
    
    double maxr, alpha, beta, sigmanew;
    
    int n = 0;
    while (n < maxIter) {
       
        A.vectorMultiply(s, z);

        alpha = sigma/(s.dot(z));
        
        p = p + alpha*s;

        r = r - alpha*z;
        
        maxr = max(r);
        
        if (maxr < tol){
            std::cout << "CG solver exiting after reaching required tolerance " << maxr << "after " << n << "iterations" << std::endl;
            return true;
        }
        
        sigmanew = r.dot(r);
        
        beta = sigmanew / sigma;
        
        s = r + beta*s;
        
        sigma = sigmanew;

        n++;

    }
    
    std::cout << "Max number of " << maxIter << "allowed iterations was exceeded, last maxDiff was " << maxr << std::endl;
    return false;
    
};


const SymmBandMatrix::VectorType operator+(const SymmBandMatrix::VectorType& lhs, const SymmBandMatrix::VectorType& rhs){
    SymmBandMatrix::VectorType result(lhs.size());
    for (int i=0; i < lhs.size(); i++){
        result[i] = lhs[i] + rhs[i];
    }
    return result;
} 


const SymmBandMatrix::VectorType operator-(const SymmBandMatrix::VectorType& lhs, const SymmBandMatrix::VectorType& rhs){
    SymmBandMatrix::VectorType result(lhs.size());
    for (int i=0; i < lhs.size(); i++){
        result[i] = lhs[i] - rhs[i];
    }
    return result;
} 


const SymmBandMatrix::VectorType operator*(const double& lhs, const SymmBandMatrix::VectorType& rhs){
    SymmBandMatrix::VectorType result{rhs};
    result *= lhs;
    return result;
}


const double max(const SymmBandMatrix::VectorType& x){
    double result = *std::max_element(x.data(), x.data() + x.size());
    return result;
}
