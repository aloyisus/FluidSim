//
//  Matrix.h
//  Fluid
//
//  Created by John Kelly on 26/01/2014.
//  Copyright (c) 2014 John Kelly. All rights reserved.
//

#ifndef __Fluid__Matrix__
#define __Fluid__Matrix__


#define matrixd matrix<double>
#define diagonald diagonal<double>
#define vectord vector<double>
#define symmbandd symm_band<double>

//#define TOL double(0.00001)
#define TOL double(0.001)

#define MAXITER 10000


#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>




template <class ElementType>
class diagonal;

template <class ElementType>
class matrix;

template <class ElementType>
class vector;

template <class ElementType>
class symm_band;


//Base class for all matrices and vectors
template <class ElementType> class base_matrix {
    
 
public:

    
    base_matrix();
    
    ~base_matrix() {
        delete[] data;
    }
    
    ElementType at(int row, int column) const;
    ElementType &at(int row, int column);

    base_matrix<ElementType> operator+(const base_matrix<ElementType>& op2);
    base_matrix<ElementType> operator*(const base_matrix<ElementType>& op2);
    base_matrix<ElementType> operator-(const base_matrix<ElementType>& op2);
    

    
    
protected:

    ElementType* data;
    
    int rows;
    int columns;
    
    
};




template <class ElementType>
class matrix : protected base_matrix<ElementType>{
    
    
public:

    matrix(int row, int column);
    
    //Initialize with value
    matrix(int row, int column, ElementType value);
        
    //Copy constructor
    matrix(const matrix &a);
    
    //Assignment operator
    matrix<ElementType>& operator=(const matrix<ElementType>& op2);
    
    void fill(ElementType f);
    ElementType at(int row, int column) const;
    ElementType& at(int row, int column);
    
    
    matrix<ElementType> transpose();
    
    
    matrix<ElementType> operator+(const matrix<ElementType>& op2);
    matrix<ElementType> operator*(const matrix<ElementType>& op2);
    matrix<ElementType> operator-(const matrix<ElementType>& op2);
    
    template <class T1>
    friend matrix<T1> operator+(const diagonal<T1>& op1, const matrix<T1>& op2);
    
    template <class T1>
    friend matrix<T1> operator+(const matrix<T1>& op1, const diagonal<T1>& op2);
    
    template <class T1>
    friend matrix<T1> operator*(const diagonal<T1>& op1, const matrix<T1>& op2);
    
    template <class T1>
    friend matrix<T1> operator*(const matrix<T1>& op1, const diagonal<T1>& op2);
    
    template <class T1>
    friend vector<T1> operator*(const matrix<T1>& op1, const vector<T1>& op2);
    
    template <class T1>
    friend vector<T1> operator*(const vector<T1>& op1, const matrix<T1>& op2);
    
    //Scalar multiplication aB = C
    template <class T1>
    friend matrix<T1> operator*(const T1 & scalar, const matrix<T1>& op2);
    
};




template <class ElementType>
class diagonal: protected base_matrix<ElementType>{
    
    
public:
    
    diagonal(int w);
    
    //Initialize with value
    diagonal(int w, ElementType value);
    
    //Copy constructor
    diagonal(const diagonal &a);
    
    //Assignment operator
    diagonal<ElementType>& operator=(const diagonal<ElementType>& op2);
    
    
    void fill(ElementType f);
    ElementType at (int r) const;
    ElementType& at(int r);
    
    diagonal<ElementType> inverse();
    
    ElementType trace();
   
    diagonal<ElementType> operator+(const diagonal<ElementType>& op2);
    diagonal<ElementType> operator*(const diagonal<ElementType>& op2);
    diagonal<ElementType> operator-(const diagonal<ElementType>& op2);
    
    template <class T1>
    friend matrix<T1> operator+(const diagonal<T1>& op1, const matrix<T1>& op2);
    
    template <class T1>
    friend matrix<T1> operator+(const matrix<T1>& op1, const diagonal<T1>& op2);
    
    template <class T1>
    friend matrix<T1> operator*(const diagonal<T1>& op1, const matrix<T1>& op2);
    
    template <class T1>
    friend matrix<T1> operator*(const matrix<T1>& op1, const diagonal<T1>& op2);
    
    template <class T1>
    friend vector<T1> operator*(const diagonal<T1>& op1, const vector<T1>& op2);
    
    template <class T1>
    friend vector<T1> operator*(const vector<T1>& op1, const diagonal<T1>& op2);
    
    //Scalar multiplication aB = C
    template <class T1>
    friend diagonal<T1> operator*(const T1 & scalar, const diagonal<T1>& op2);
    
    template <class T1>
    friend void applyPreconditioner(const symm_band<T1>& A, const diagonal<T1>& Ei, const vector<T1>& r, vector<T1>& z);
    
    template <class T1>
    friend void PCGSolver(const symm_band<T1>& A, const diagonal<T1>& Ei, const vector<T1>& b, vector<T1>& p, T1 tol, int maxIter);
    
    template <class T1>
    friend void MIC0precon(const symm_band<T1>& A, diagonal<T1>& Ei);
    
   
};




template <class ElementType>
class vector: protected base_matrix<ElementType>{
    
    
public:
    
    vector(int w);
    
    //Initialize with value
    vector(int w, ElementType value);
    
    //Copy constructor
    vector(const vector &a);
    
    //Assignment operator
    vector<ElementType>& operator=(const vector<ElementType>& op2);
    
    
    void fill(ElementType f);
    ElementType at(int row) const;
    ElementType &at(int row);
   
    ElementType max();
    ElementType sum();
    
    int getRows() {  return this->rows;  }
    
    
    vector<ElementType> operator+(const vector<ElementType>& op2);
    
    //Multiplication of two vectors is dot product
    ElementType operator*(const vector<ElementType>& op2);
    
    vector<ElementType> operator-(const vector<ElementType>& op2) const;
    
    template <class T1>
    friend vector<T1> operator*(const diagonal<T1>& op1, const vector<T1>& op2);
    
    template <class T1>
    friend vector<T1> operator*(const vector<T1>& op1, const diagonal<T1>& op2);
    
    template <class T1>
    friend vector<T1> operator*(const matrix<T1>& op1, const vector<T1>& op2);
    
    template <class T1>
    friend vector<T1> operator*(const vector<T1>& op1, const matrix<T1>& op2);
    
    //Scalar multiplication aB = C
    template <class T1>
    friend vector<T1> operator*(const T1 & scalar, const vector<T1>& op2);
    
    template <class T1>
    friend vector<T1> operator*(const symm_band<T1>& op1, const vector<T1>& op2);
    
    //Linear solver
    template <class T1>
    friend void GaussSeidelSolver(const symm_band<T1>& op1, const vector<T1>& op2, vector<T1>& b, T1 tol, int maxIter);
    
    template <class T1>
    friend void PCGSolver(const symm_band<T1>& A, const diagonal<T1>& Ei, const vector<T1>& b, vector<T1>& p, T1 tol, int maxIter);
    
    template <class T1>
    friend void CGSolver(const symm_band<T1>& A, const vector<T1>& b, vector<T1>& p, T1 tol, int maxIter);
    
    template <class T1>
    friend void applyPreconditioner(const symm_band<T1>& A, const diagonal<T1>& Ei, const vector<T1>& r, vector<T1>& z);
    
    
};




//Sparse matrix class - specifically Symmetric Band matrix, used for Bridson's A matrix and operations
//
//This is suitable for a symmetric sparse matrix with non-zero leading diag
template <class ElementType>
class symm_band : protected base_matrix<ElementType>{
    
    
public:
    
    symm_band(int size, int bandwidth, int offset1, int offset2);
    symm_band(int size, int bandwidth, int offset1, int offset2, int offset3);
    
    //Initialize with value
    symm_band(int size, int bandwidth, int offset1, int offset2, ElementType value);
    symm_band(int size, int bandwidth, int offset1, int offset2, int offset3, ElementType value);
    
    //Copy constructor
    symm_band(const symm_band &a);
    
    //Assignment operator =
    symm_band<ElementType>& operator=(const symm_band<ElementType>& op2);
 
    //Comparison equality operator ==
    symm_band<ElementType>& operator==(const symm_band<ElementType>& op2);
    
    void fill(ElementType value);
    
    ElementType max();
    
    bool testDiagonal();
    
    ElementType at(int diag, int offset) const;
    ElementType &at(int diag, int offset);
    
    //Scalar multiplication aB = C
    template <class T1>
    friend symm_band<T1> operator*(const T1 & scalar, const symm_band<T1>& op2);
    
    template <class T1>
    friend vector<T1> operator*(const symm_band<T1>& op1, const vector<T1>& op2);
    
    //Linear Solver
    template <class T1>
    friend void GaussSeidelSolver(const symm_band<T1>& op1, const vector<T1>& op2, vector<T1>& b, T1 tol, int maxIter);
    
    template <class T1>
    friend void PCGSolver(const symm_band<T1>& A, const diagonal<T1>& Ei, const vector<T1>& b, vector<T1>& p, T1 tol, int maxIter);
    
    template <class T1>
    friend void CGSolver(const symm_band<T1>& A, const vector<T1>& b, vector<T1>& p, T1 tol, int maxIter);
    
    template <class T1>
    friend void applyPreconditioner(const symm_band<T1>& A, const diagonal<T1>& Ei, const vector<T1>& r, vector<T1>& z);
    
    template <class T1>
    friend void MIC0precon(const symm_band<T1>& A, diagonal<T1>& Ei);
    

protected:
    
    int bandwidth;
    int offset1;
    int offset2;
    int offset3;

};





#endif /* defined(__Fluid__Matrix__) */
