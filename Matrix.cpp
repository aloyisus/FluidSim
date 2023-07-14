//
//  Matrix.cpp
//  Fluid
//
//  Created by John Kelly on 26/01/2014.
//  Copyright (c) 2014 John Kelly. All rights reserved.
//

#include "matrix.h"




/*
//
//    base_matrix definitions
//
//
*/

template <class ElementType>
base_matrix<ElementType>::base_matrix(){
    
    
};




/*
 //
 //    matrix definitions
 //
 //
 */

template <class ElementType>
matrix<ElementType>::matrix(int row, int column){
    
    this->rows = row;
    this->columns = column;
    
    this->data = new ElementType [row*column];
    
    
};


//Initialize with value
template <class ElementType>
matrix<ElementType>::matrix(int row, int column, ElementType value){
    
    this->rows = row;
    this->columns = column;
    
    this->data = new ElementType [row*column];

    for(int i=0; i<row*column; i++)
        this->data[i] = value;
    

    
};


//Copy constructor
template <class ElementType>
matrix<ElementType>::matrix(const matrix &a){
    
    this->rows = a.rows;
    this->columns = a.columns;
    
    this->data = new ElementType [a.rows*a.columns];
    
    for(int i=0; i<a.rows*a.columns; i++)
        this->data[i] = a.data[i];
    
    
};


//Assignment operator
template <class ElementType>
matrix<ElementType>& matrix<ElementType>::operator=(const matrix<ElementType>& op2){
    
    // check for self-assignment
    if (this == &op2)
        return *this;
    
    this->rows=op2.rows;
    this->columns = op2.columns;
    
    for(int i=0; i<op2.rows*op2.columns; i++){
        
        this->data[i] = op2.data[i];
       
    }
    
    return *this;
    
};


template <class ElementType>
ElementType matrix<ElementType>::at (int r,int c) const {
    
    //elements with the same row number are stored contiguously in memory
    return this->data[this->columns*r+c];
    
};


template <class ElementType>
ElementType& matrix<ElementType>::at(int r,int c){
    
    //elements with the same row number are stored contiguously in memory
    return this->data[this->columns*r+c];
    
};


template <class ElementType>
void matrix<ElementType>::fill(ElementType f)
{
    
    for(int i=0; i<this->rows*this->columns; i++){
        
        this->data[i] = f;
        
    }
    
};


template <class ElementType>
matrix<ElementType> matrix<ElementType>::transpose(){

    int row = this->rows;
    int column = this->columns;
    
    matrix<ElementType> temp = matrix<ElementType>(row,column);
    
    //A.transpose(rc) = A(cr)
    //
    //r = i/column, c = i%column
    //
    //A(cr) = A[(c*column + r]
    for(int i=0; i<row*column; i++){
        
        temp.data[i] = this->data[(i%column)*column+(i/column)];
        
    }
    
    
    return temp;
    
};


template <class ElementType>
matrix<ElementType> matrix<ElementType>::operator+(const matrix<ElementType> &op2){
    
    matrix<ElementType> temp = matrix<ElementType>(this->rows,this->columns);
        
    
    for(int i=0; i<this->rows*this->columns; i++){
 
        temp.data[i] = this->data[i] + op2.data[i];
       
    }
    
    return temp;
    
}


template <class ElementType>
matrix<ElementType> matrix<ElementType>::operator*(const matrix<ElementType>& op2){
    
    int out_columns = op2.columns;
    int out_rows = this->rows;
    matrix<ElementType> temp = matrix<ElementType>(out_rows,out_columns,0);
    
    for(int r=0; r<out_columns; r++)
        for(int c=0; c<out_rows; c++){
            
            //A(rc) = sum over i:  B(ri)C(ic)
            for(int i=0; i<this->columns; i++){
                
                temp.data[r*out_columns+c] += this->data[r*this->columns+i] * op2.data[i*out_columns+c];
                
            }
            
        }
    
    return temp;
    
}


template <class ElementType>
matrix<ElementType> matrix<ElementType>::operator-(const matrix<ElementType>& op2){
    
    matrix<ElementType> temp = matrix<ElementType>(this->rows,this->columns);
    
    
    for(int i=0; i<this->rows*this->columns; i++){
        
        temp.data[i] = this->data[i] - op2.data[i];
        
    }
    
    return temp;
    
}











/*
 //
 //    diagonal definitions
 //
 //
 */

template <class ElementType>
diagonal<ElementType>::diagonal(int w){
    
    diagonal::rows = w;
    diagonal::columns = w;
    
    diagonal::data = new ElementType [w];
    
    
};


//Initialize with value
template <class ElementType>
diagonal<ElementType>::diagonal(int w, ElementType value){
    
    this->rows = w;
    this->columns = w;
    
    this->data = new ElementType [w];
    
    for(int i=0; i<w; i++)
        this->data[i] = value;
    
    
    
};


//Copy constructor
template <class ElementType>
diagonal<ElementType>::diagonal(const diagonal &a){
    
    
    this->rows = a.rows;
    this->columns = a.columns;
    
    this->data = new ElementType [a.rows];
    
    for(int i=0; i<a.rows; i++)
        this->data[i] = a.data[i];
    
    
};


//Assignment operator
template <class ElementType>
diagonal<ElementType>& diagonal<ElementType>::operator=(const diagonal<ElementType>& op2){
    
    // check for self-assignment
    if (this == &op2)
        return *this;
    
    this->rows=op2.rows;
    this->columns = op2.columns;
    
    for(int i=0; i<op2.rows; i++){
        
        this->data[i] = op2.data[i];
        
    }
    
    return *this;
    
};


template <class ElementType>
ElementType diagonal<ElementType>::at (int r) const {
    
    //elements with the same row number are stored contiguously in memory
    return this->data[r];
    
};


template <class ElementType>
ElementType& diagonal<ElementType>::at(int r){
    
    //elements with the same row number are stored contiguously in memory
    return this->data[r];
    
};


template <class ElementType>
void diagonal<ElementType>::fill(ElementType f){
    
    for(int i=0; i<diagonal::rows; i++){
        
        diagonal::data[i] = f;
        
    }
    
    
};


template <class ElementType>
diagonal<ElementType> diagonal<ElementType>::inverse(){

    int row = this->rows;
    diagonal<ElementType> temp = diagonal<ElementType>(row);
    
    //A.inverse(ij) = 1/A(ij)
    for(int i=0; i<row; i++){
        
        temp.data[i] = 1 / this->data[i];
        
    }
    
    return temp;
    
};


template <class ElementType>
ElementType diagonal<ElementType>::trace(){
    
    ElementType temp = 0;
    
    //A.trace = sum over i:  A(ii)     this is only defined for square matrices, like this diagonal matrix
    for(int i=0; i<this->rows; i++){
        
        //int d =i*this->columns + i;
        temp += this->data[i];
        
    }
    
    return temp;
    
};


template <class ElementType>
diagonal<ElementType> diagonal<ElementType>::operator+(const diagonal<ElementType>& op2){
    
    int row = this->rows;
    diagonal<ElementType> temp = diagonal<ElementType>(row);
    
    
    for(int i=0; i<row; i++){
        
        temp.data[i] = this->data[i] + op2.data[i];
        
    }
    
    
    return temp;
    
};


template <class ElementType>
diagonal<ElementType> diagonal<ElementType>::operator*(const diagonal<ElementType>& op2){
    
    int row = this->rows;
    diagonal<ElementType> temp = diagonal<ElementType>(row);
    
    
    //A(rr) = sum over i:  B(ri)C(ir) = B(rr)C(rr)     (using the fact A is diagonal)
    for(int i=0; i<row; i++){
        
        temp.data[i] = this->data[i] * op2.data[i];
        
    }
    
    
    return temp;
    
};


template <class ElementType>
diagonal<ElementType> diagonal<ElementType>::operator-(const diagonal<ElementType>& op2){
    
    int row = this->rows;
    diagonal<ElementType> temp = diagonal<ElementType>(row);
    
    
    for(int i=0; i<row; i++){
        
        temp.data[i] = this->data[i] - op2.data[i];
        
    }
    
    
    return temp;
    
};






/*
 //
 //    vector definitions
 //
 //
 */

template <class ElementType>
vector<ElementType>::vector(int w){
    
    vector::rows = w;
    vector::columns = 1;
    
    vector::data = new ElementType [w];
    
    
};


//Initialize with value
template <class ElementType>
vector<ElementType>::vector(int w, ElementType value){
    
    vector::rows = w;
    vector::columns = 1;
    
    vector::data = new ElementType [w];
    
    for(int i=0; i<w; i++)
        vector::data[i] = value;
    
    
    
};


//Copy constructor
template <class ElementType>
vector<ElementType>::vector(const vector &a){
    
    
    int w = a.rows;
    
    this->rows = w;
    this->columns = a.columns;
    
    this->data = new ElementType [w];
    
    for(int i=0; i<w; i++)
        this->data[i] = a.data[i];
    
    
};


//Assignment operator
template <class ElementType>
vector<ElementType>& vector<ElementType>::operator=(const vector<ElementType>& op2){
    
    // check for self-assignment
    if (this == &op2)
        return *this;
    
    this->rows=op2.rows;
    this->columns = op2.columns;
    
    for(int i=0; i<op2.rows; i++){
        
        this->data[i] = op2.data[i];
        
    }
    
    return *this;
    
};


template <class ElementType>
void vector<ElementType>::fill(ElementType f){
    
    for(int i=0; i<vector::rows; i++){
        
        vector::data[i] = f;
        
    }
    
    
};


template <class ElementType>
ElementType vector<ElementType>::max(){
    
    int size = this->rows;
    ElementType max = 0;
    
    for(int i=0; i<size; i++)
        max = fmax(max,fabs(this->data[i]));
    
    return max;
    
};


template <class ElementType>
ElementType vector<ElementType>::sum(){
    
    int size = this->rows;
    ElementType sum = 0;
    
    for(int i=0; i<size; i++)
        sum += this->data[i];
    
    return sum;
    
};


template <class ElementType>
ElementType vector<ElementType>::at (int r) const {
    
    //elements with the same row number are stored contiguously in memory
    return this->data[r];
    
};


template <class ElementType>
ElementType& vector<ElementType>::at(int r){
    
    //elements with the same row number are stored contiguously in memory
    return this->data[r];
    
};


//Multiplication of two vectors is dot product
template <class ElementType>
ElementType vector<ElementType>::operator*(const vector<ElementType>& op2){
    
    int w = op2.rows;
    ElementType temp = 0;
    
    
    for(int i=0; i<w; i++){
        
        temp += this->data[i] * op2.data[i];
        
    }
    
    
    return temp;
    
};


template <class ElementType>
vector<ElementType> vector<ElementType>::operator+(const vector<ElementType>& op2){
    
    int w = op2.rows;
    vector<ElementType> temp = vector<ElementType>(w);
    
    
    for(int i=0; i<w; i++){
        
        temp.data[i] = this->data[i] + op2.data[i];
        
    }
    
    
    return temp;
    
};


template <class ElementType>
vector<ElementType> vector<ElementType>::operator-(const vector<ElementType>& op2) const{
    
    int w = op2.rows;
    vector<ElementType> temp = vector<ElementType>(w);
    
    
    for(int i=0; i<w; i++){
        
        temp.data[i] = this->data[i] - op2.data[i];
        
    }
    
    
    return temp;
    
};








/*
 //
 //    symm_band definitions
 //
 //
 */

template <class ElementType>
symm_band<ElementType>::symm_band(int size, int bandwidth, int offset1, int offset2){
    
    symm_band::rows = size;
    
    symm_band::bandwidth = bandwidth;
    symm_band::offset1 = offset1;
    symm_band::offset2 = offset2;
    
    int column = (bandwidth+1)/2;
    
    symm_band::columns = column;
    
    symm_band::data = new ElementType[size*column];
    
    
};

template <class ElementType>
symm_band<ElementType>::symm_band(int size, int bandwidth, int offset1, int offset2, int offset3){
    
    symm_band::rows = size;
    
    symm_band::bandwidth = bandwidth;
    symm_band::offset1 = offset1;
    symm_band::offset2 = offset2;
    symm_band::offset3 = offset3;
    
    int column = (bandwidth+1)/2;
    
    symm_band::columns = column;
    
    symm_band::data = new ElementType[size*column];
    
    
};


//Copy constructor
template <class ElementType>
symm_band<ElementType>::symm_band(const symm_band &a){
    
    
    int size = a.rows;
    this->rows = size;
    
    int column = a.columns;
    this->columns = column;
    
    int bandwidth = a.bandwidth;
    this->bandwidth = bandwidth;
    
    this->offset1 = a.offset1;
    this->offset2 = a.offset2;
    this->offset3 = a.offset3;
    
    int arraylength=size*column;
    symm_band::data = new ElementType [arraylength];
    
    for(int i=0; i<arraylength; i++)
        this->data[i] = a.data[i];
    
};


template <class ElementType>
symm_band<ElementType>& symm_band<ElementType>::operator==(const symm_band<ElementType>& op2){
    
    if (this->rows != op2.rows)
        return false;
    if (this->columns != op2.columns)
        return false;
    if (this->bandwidth != op2.bandwidth)
        return false;
    if (this->offset1 != op2.offset1)
        return false;
    if (this->offset2 != op2.offset2)
        return false;
    if (this->offset3 != op2.offset3)
        return false;
        
    int row = this->rows;
    int column = this->columns;
    
    for(int i=0; i<row*column; i++)
        if (this->data[i] != op2.data[i])
            return false;
    
    return true;
    
}


template <class ElementType>
void symm_band<ElementType>::fill(ElementType value){
    
    
    int rows = this->rows;
    int columns = this->columns;
    
    int arraylength=rows*columns;
    
    for(int i=0; i<arraylength; i++){
        
        this->data[i] = value;
        
    }
   
}

template <class ElementType>
ElementType symm_band<ElementType>::max(){
    
    ElementType tempmax = 0;
    
    int rows = this->rows;
    int columns = this->columns;
    
    int arraylength=rows*columns;
    
    for(int i=0; i<arraylength; i++){
        
        tempmax = std::max(fabs(this->data[i]),tempmax);
        
    }
   
    return tempmax;
}


template <class ElementType>
bool symm_band<ElementType>::testDiagonal(){
    
    int rows = this->rows;
    int columns = this->columns;
    
    for(int i=0; i<rows; i++){
        
        if (this->data[i*columns]==0)
            return false;
        
    }
    
    return true;
    
}


//Assignment operator
template <class ElementType>
symm_band<ElementType>& symm_band<ElementType>::operator=(const symm_band<ElementType>& op2){
    
    // check for self-assignment
    if (this == &op2)
        return *this;
    
    int size = op2.rows;
    this->rows = size;
    
    int column = op2.columns;
    this->columns = column;
    
    int bandwidth = op2.bandwidth;
    this->bandwidth = bandwidth;
    
    this->offset1 = op2.offset1;
    this->offset2 = op2.offset2;
    this->offset3 = op2.offset3;
    
    int arraylength=size*column;
    
    for(int i=0; i<arraylength; i++)
        this->data[i] = op2.data[i];
    
    return *this;
    
};


template <class ElementType>
ElementType& symm_band<ElementType>::at(int diag, int offset) {
    
    
    //  This is how the matrix is arranged. It is symmetric, so we
    //  don't need to keep a little less than half the values, i.e. a32 = a23
    //
    //  Notice the first column contains the diagonal elements always
    //
    //  The width of the stored data is symm_band::columns,
    //  this is not the same as the width of the matrix we are
    //  representing, which is a square matrix of symm_band::rows * symm_band::rows
    //
    //  This code does no checking to see if the offset (or the diagonal index) is out of bounds
    //
    //
    //      [a00,   a01,   a02,   a03]
    //      [a11,   a12,   a13,   a24]
    //      [a22,   a23,   a24,   a25]
    //      [a33,   a34,   a35,   a36]
    //      [...                     ]
    //      [aii, aii+1, aii+2, aii+3]
    //      [                        ]
    //
    //
    //
    //
        
    return this->data[diag*this->columns+offset];
    
};


template <class ElementType>
ElementType symm_band<ElementType>::at(int diag, int offset) const {
    
    
    //  This is how the matrix is arranged. It is symmetric, so we
    //  don't need to keep a little less than half the values, i.e. a32 = a23
    //
    //  Notice the first column contains the diagonal elements always
    //
    //  The width of the stored data is symm_band::columns,
    //  this is not the same as the width of the matrix we are
    //  representing, which is a square matrix of symm_band::rows * symm_band::rows
    //
    //  This code does no checking to see if the offset (or the diagonal index) is out of bounds
    //
    //
    //      [a00,   a01,   a02,   a03]
    //      [a11,   a12,   a13,   a24]
    //      [a22,   a23,   a24,   a25]
    //      [a33,   a34,   a35,   a36]
    //      [...                     ]
    //      [aii, aii+1, aii+2, aii+3]
    //      [                        ]
    //
    //
    //
    //
    
    return this->data[diag*this->columns+offset];
    
};



    








//
//
//FRIEND FUNCTIONS FOR BINARY OPERATIONS BETWEEN DIFFERENT CLASSES OF MATRIX
//
//
template <class T1>
matrix<T1> operator+(diagonal<T1>& op1, matrix<T1>& op2){
    
    
    matrix<T1> temp = op2;
    
    
    for(int i=0; i<op1.rows; i++){
        
        temp.data[i*temp.columns+i] += op1.data[i];
        
    }
    
    
    return temp;
    
    
};


template <class T1>
matrix<T1> operator+(const matrix<T1>& op1, const diagonal<T1>& op2){
    
    
    matrix<T1> temp = op1;
    
    
    for(int i=0; i<op2.rows; i++){
        
        temp.data[i*temp.columns+i] += op2.data[i];
        
    }
    
    
    return temp;
    
    
};


template <class T1>
matrix<T1> operator*(diagonal<T1>& op1, matrix<T1>& op2){
    
    int row2 = op2.rows;
    int column2 = op2.columns;
    
    //Diagonal matrix (always square) must have rows=columns = no. of rows in op2, i.e. row2 x row2
    //there is no code here to check this
    
    //Result matrix will be same size as op2, i.e. row2 * column2
    matrix<T1> temp = op2;
    
    //A(rc) = B(rr)C(rc)     where B is the diagonal matrix, and result A is not necessarily diagonal
    //
    //B(rr) = B[i/column2]
    for(int i=0; i<row2*column2; i++){
        
        temp.data[i] *= op1.data[i/column2];
        
    }
    
    return temp;
    
    
};


template <class T1>
matrix<T1> operator*(const matrix<T1>& op1, const diagonal<T1>& op2){
    
    int row1 = op1.rows;
    int column1 = op1.columns;
    
    //Diagonal matrix (always square) must have rows=columns = no. of columns in op1, i.e. column1 x column1
    //there is no code here to check this
    
    //Result matrix will be same size as op1, i.e. row1 * column1
    matrix<T1> temp = op1;
    
    
    //A(rc) = B(cc)C(rc)     where B is the diagonal matrix, and result A is not necessarily diagonal
    //
    //B(cc) = B[i%column1]      % gives the remainder
    for(int i=0; i<row1*column1; i++){
        
        temp.data[i] *= op2.data[i%column1];
        
    }
    
    return temp;
    
    
};


//Scalar multiplication aB = C
template <class T1>
matrix<T1> operator*(const T1 & scalar, const matrix<T1>& op2){
    
    matrix<T1> temp = op2;
    
    for(int i=0; i<op2.rows*op2.columns; i++){
        
        temp.data[i] *= scalar;
        
    }
    
    return temp;
    
    
};


//Scalar multiplication aB = C
template <class T1>
diagonal<T1> operator*(const T1 & scalar, const diagonal<T1>& op2){
    
    diagonal<T1> temp = op2;
    
    for(int i=0; i<op2.rows; i++){
        
        temp.data[i] *= scalar;
        
    }
    
    return temp;
    
    
};


//Scalar multiplication aB = C
template <class T1>
vector<T1> operator*(const T1 & scalar, const vector<T1>& op2){
    
    vector<T1> temp = op2;
    
    for(int i=0; i<op2.rows; i++){
        
        temp.data[i] *= scalar;
        
    }
    
    return temp;
    
    
};


template <class T1>
vector<T1> operator*(const diagonal<T1>& op1, const vector<T1>& op2){

    int column1 = op1.columns;
    
    //Diagonal matrix (always square) must have rows=columns = no. of columns in op1, i.e. column1 x column1
    //Vector op2 must be a column matrix with rows = number of columns in op1 = column1
    //there is no code here to check either of these things
    
    //Result vector will be same size as op2
    vector<T1> temp = op2;

    
    for(int i=0; i<column1; i++){
        
        temp.data[i] *= op1.data[i];
        
    }
    
    return temp;
    
};


template <class T1>
vector<T1> operator*(const vector<T1>& op1, const diagonal<T1>& op2){
    
    int column1 = op2.columns;
    
    //Diagonal matrix (always square) must have rows=columns = no. of columns in op2, i.e. column1 x column1
    //Vector op1 must be a column matrix with rows = number of columns in op1 = column1
    //(actually, it should be a row matrix, but the maths is unaffected by treating row and column matrices as the same thing,
    //what I'm basically doing is, for column matrix a and row matrix a^T, treating a^T = a in the code.
    //There is no code here to check if the operands are of the correct shape
    
    //Note we get exactly the same result as Dv i.e. v^T*D = D*v. To see this, (D*v)^T = v^T*D^T = v^T*D
    
    //Result vector will be same size as op1
    vector<T1> temp = op1;
    
    
    for(int i=0; i<column1; i++){
        
        temp.data[i] *= op1.data[i];
        
    }
    
    return temp;
    
};


template <class T1>
vector<T1> operator*(const matrix<T1>& op1, const vector<T1>& op2){
    
    //a = Ab (A=op1, b=op2)
    //b must have the same number of rows as A has columns, a will have same number of rows as A, and 1 column
    
    int row1 = op1.rows;
    int column1 = op1.columns;

    vector<T1> temp = vector<T1>(row1, 0);
    
    for(int r=0; r<row1; r++){
        for(int i=0; i<column1; i++){
            
            temp.data[r] += op1.data[r*column1+i]*op2.data[i];
        
        }
    }
    
    return temp;
    
};


template <class T1>
vector<T1> operator*(const vector<T1>& op1, const matrix<T1>& op2){
    
    //a^T = b^TA (A=op2, b^T=op1)
    //b^T must have the same number of columns as A has rows, a^T will have same number of columns as A, and 1 row
    
    //Further note - a and a^T are treated identically by the code - both are column vectors
    
    int row2 = op2.rows;
    int column2 = op2.columns;
    
    vector<T1> temp = vector<T1>(column2, 0);
    
    for(int c=0; c<column2; c++){
        for(int i=0; i<row2; i++){
            
            temp.data[c] += op1.data[i]*op2.data[i*column2+c];
            
        }
    }
    
    return temp;
    
};


template <class T1>
vector<T1> operator*(const symm_band<T1>& op1, const vector<T1>& op2){

	int size = op1.rows;
	int column = op1.columns;
    
    int off1 = op1.offset1;
    int off2 = op1.offset2;
    int off3 = op1.offset3;
   
    
    vector<T1> temp(size,0);
	
    
    switch(op1.bandwidth){
            
        case 5:

            for(int i=0; i<off1; i++){
                
                //diagonal
                temp.data[i] += op1.data[column*i] * op2.data[i];
                
                //offset1 to right of diagonal
                temp.data[i] += op1.data[column*i+1] * op2.data[i+off1];
                
                //offset2 to right of diagonal
                temp.data[i] += op1.data[column*i+2] * op2.data[i+off2];
                
            }
            
            for(int i=off1; i<off2; i++){
                
                //offset1 to left of diagonal
                temp.data[i] += op1.data[column*(i-off1)+1] * op2.data[i-off1];
                
                //diagonal
                temp.data[i] += op1.data[column*i] * op2.data[i];
                
                //offset1 to right of diagonal
                temp.data[i] += op1.data[column*i+1] * op2.data[i+off1];
                
                //offset2 to right of diagonal
                temp.data[i] += op1.data[column*i+2] * op2.data[i+off2];
                
            }
            
            for(int i=off2; i<size-off2; i++){
                
                //offset2 to left of diagonal
                temp.data[i] += op1.data[column*(i-off2)+2] * op2.data[i-off2];
                
                //offset1 to left of diagonal
                temp.data[i] += op1.data[column*(i-off1)+1] * op2.data[i-off1];
                
                //diagonal
                temp.data[i] += op1.data[column*i] * op2.data[i];
                
                //offset1 to right of diagonal
                temp.data[i] += op1.data[column*i+1] * op2.data[i+off1];
                
                //offset2 to right of diagonal
                temp.data[i] += op1.data[column*i+2] * op2.data[i+off2];
                
            }
            
            for(int i=size-off2; i<size-off1; i++){
                
                //offset2 to left of diagonal
                temp.data[i] += op1.data[column*(i-off2)+2] * op2.data[i-off2];
                
                //offset1 to left of diagonal
                temp.data[i] += op1.data[column*(i-off1)+1] * op2.data[i-off1];
                
                //diagonal
                temp.data[i] += op1.data[column*i] * op2.data[i];
                
                //offset1 to right of diagonal
                temp.data[i] += op1.data[column*i+1] * op2.data[i+off1];
                
            }
            
            for(int i=size-off1; i<size; i++){
                
                //offset2 to left of diagonal
                temp.data[i] += op1.data[column*(i-off2)+2] * op2.data[i-off2];
                
                //offset1 to left of diagonal
                temp.data[i] += op1.data[column*(i-off1)+1] * op2.data[i-off1];
                
                //diagonal
                temp.data[i] += op1.data[column*i] * op2.data[i];
                
            }
    
            break;
            
        case 7:
            

            for(int i=0; i<off1; i++){
                
                //diagonal
                temp.data[i] += op1.data[column*i] * op2.data[i];
                
                //offset1 to right of diagonal
                temp.data[i] += op1.data[column*i+1] * op2.data[i+off1];
                
                //offset2 to right of diagonal
                temp.data[i] += op1.data[column*i+2] * op2.data[i+off2];
                
                //offset3 to right of diagonal
                temp.data[i] += op1.data[column*i+3] * op2.data[i+off3];
                
            }
            
            for(int i=off1; i<off2; i++){
                
                //offset1 to left of diagonal
                temp.data[i] += op1.data[column*(i-off1)+1] * op2.data[i-off1];
                
                //diagonal
                temp.data[i] += op1.data[column*i] * op2.data[i];
                
                //offset1 to right of diagonal
                temp.data[i] += op1.data[column*i+1] * op2.data[i+off1];
                
                //offset2 to right of diagonal
                temp.data[i] += op1.data[column*i+2] * op2.data[i+off2];
                
                //offset3 to right of diagonal
                temp.data[i] += op1.data[column*i+3] * op2.data[i+off3];
                
            }
            
            for(int i=off2; i<off3; i++){
                
                //offset2 to left of diagonal
                temp.data[i] += op1.data[column*(i-off2)+2] * op2.data[i-off2];
                
                //offset1 to left of diagonal
                temp.data[i] += op1.data[column*(i-off1)+1] * op2.data[i-off1];
                
                //diagonal
                temp.data[i] += op1.data[column*i] * op2.data[i];
                
                //offset1 to right of diagonal
                temp.data[i] += op1.data[column*i+1] * op2.data[i+off1];
                
                //offset2 to right of diagonal
                temp.data[i] += op1.data[column*i+2] * op2.data[i+off2];
                
                //offset3 to right of diagonal
                temp.data[i] += op1.data[column*i+3] * op2.data[i+off3];
                
            }
            
            for(int i=off3; i<size-off3; i++){
                
                //offset3 to left of diagonal
                temp.data[i] += op1.data[column*(i-off3)+3] * op2.data[i-off3];
                
                //offset2 to left of diagonal
                temp.data[i] += op1.data[column*(i-off2)+2] * op2.data[i-off2];
                
                //offset1 to left of diagonal
                temp.data[i] += op1.data[column*(i-off1)+1] * op2.data[i-off1];
                
                //diagonal
                temp.data[i] += op1.data[column*i] * op2.data[i];
                
                //offset1 to right of diagonal
                temp.data[i] += op1.data[column*i+1] * op2.data[i+off1];
                
                //offset2 to right of diagonal
                temp.data[i] += op1.data[column*i+2] * op2.data[i+off2];
                
                //offset3 to right of diagonal
                temp.data[i] += op1.data[column*i+3] * op2.data[i+off3];
                
            }
            
            for(int i=size-off3; i<size-off2; i++){
                
                //offset3 to left of diagonal
                temp.data[i] += op1.data[column*(i-off3)+3] * op2.data[i-off3];
                
                //offset2 to left of diagonal
                temp.data[i] += op1.data[column*(i-off2)+2] * op2.data[i-off2];
                
                //offset1 to left of diagonal
                temp.data[i] += op1.data[column*(i-off1)+1] * op2.data[i-off1];
                
                //diagonal
                temp.data[i] += op1.data[column*i] * op2.data[i];
                
                //offset1 to right of diagonal
                temp.data[i] += op1.data[column*i+1] * op2.data[i+off1];
                
                //offset2 to right of diagonal
                temp.data[i] += op1.data[column*i+2] * op2.data[i+off2];
                
            }
            
            for(int i=size-off2; i<size-off1; i++){
                
                //offset3 to left of diagonal
                temp.data[i] += op1.data[column*(i-off3)+3] * op2.data[i-off3];
                
                //offset2 to left of diagonal
                temp.data[i] += op1.data[column*(i-off2)+2] * op2.data[i-off2];
                
                //offset1 to left of diagonal
                temp.data[i] += op1.data[column*(i-off1)+1] * op2.data[i-off1];
                
                //diagonal
                temp.data[i] += op1.data[column*i] * op2.data[i];
                
                //offset1 to right of diagonal
                temp.data[i] += op1.data[column*i+1] * op2.data[i+off1];
                
            }
            
            for(int i=size-off1; i<size; i++){
                
                //offset3 to left of diagonal
                temp.data[i] += op1.data[column*(i-off3)+3] * op2.data[i-off3];
                
                //offset2 to left of diagonal
                temp.data[i] += op1.data[column*(i-off2)+2] * op2.data[i-off2];
                
                //offset1 to left of diagonal
                temp.data[i] += op1.data[column*(i-off1)+1] * op2.data[i-off1];
                
                //diagonal
                temp.data[i] += op1.data[column*i] * op2.data[i];
                
            }
            
            break;
            
    }
    
    return temp;
    
};


//Linear Solver
template <class T1>
void GaussSeidelSolver(const symm_band<T1>& op1, const vector<T1>& op2, vector<T1>& b, T1 tol, int maxIter){
    
    //The equation we are solving is op1 * b = op2 where op1 is a symmetric band matrix, op2 is a vector,
    //and b is the vector we are trying to find
    
    int size = op1.rows;
	int column = op1.columns;
    
    int off1 = op1.offset1;
    int off2 = op1.offset2;
    int off3 = op1.offset3;
    
    T1 offDiag;
    T1 maxDiff;
    T1 oldb;
    
    //keep track of iteration number
    int n = 0;
    
    switch(op1.bandwidth){
            
        case 5:
            
            //loop until error is acceptable or we have exceeded max. iterations
            while (n < maxIter) {
                
                maxDiff = 0;
                for(int i=0; i<off1; i++){
                    
                    offDiag = 0;
                    oldb = b.data[i];
                    
                    //offset1 to right of diagonal
                    offDiag += op1.data[column*i+1] * b.data[i+off1];
                    
                    //offset2 to right of diagonal
                    offDiag += op1.data[column*i+2] * b.data[i+off2];
                    
                    //b(i) = (op2(i) - offDiag(i))/diag(i)
                    //
                    //where offDiag is all the offdiagonal elements of op1 multiplied by corresponding elements of b, diag(i) is leading diagonal element of op1
                    b.data[i] = (op2.data[i] - offDiag)/op1.data[column*i];
                    maxDiff = std::max(maxDiff, double(fabs(b.data[i] - oldb)));
                   
                }
                
                for(int i=off1; i<off2; i++){
                    
                    offDiag = 0;
                    oldb = b.data[i];
                    
                    //offset1 to left of diagonal
                    offDiag += op1.data[column*(i-off1)+1] * b.data[i-off1];
                    
                    //offset1 to right of diagonal
                    offDiag += op1.data[column*i+1] * b.data[i+off1];
                    
                    //offset2 to right of diagonal
                    offDiag += op1.data[column*i+2] * b.data[i+off2];
                    
                    //b(i) = (op2(i) - offDiag(i))/diag(i)
                    //
                    //where offDiag is all the offdiagonal elements of op1 multiplied by corresponding elements of b, diag(i) is leading diagonal element of op1
                    b.data[i] = (op2.data[i] - offDiag)/op1.data[column*i];
                    maxDiff = std::max(maxDiff, double(fabs(b.data[i] - oldb)));
                    
                }
                
                for(int i=off2; i<size-off2; i++){
                    
                    offDiag = 0;
                    oldb = b.data[i];
                    
                    //offset2 to left of diagonal
                    offDiag += op1.data[column*(i-off2)+2] * b.data[i-off2];
                    
                    //offset1 to left of diagonal
                    offDiag += op1.data[column*(i-off1)+1] * b.data[i-off1];
                    
                    //offset1 to right of diagonal
                    offDiag += op1.data[column*i+1] * b.data[i+off1];
                    
                    //offset2 to right of diagonal
                    offDiag += op1.data[column*i+2] * b.data[i+off2];
                    
                    //b(i) = (op2(i) - offDiag(i))/diag(i)
                    //
                    //where offDiag is all the offdiagonal elements of op1 multiplied by corresponding elements of b, diag(i) is leading diagonal element of op1
                    b.data[i] = (op2.data[i] - offDiag)/op1.data[column*i];
                    maxDiff = std::max(maxDiff, double(fabs(b.data[i] - oldb)));
                    
                }
                
                for(int i=size-off2; i<size-off1; i++){
                    
                    offDiag = 0;
                    oldb = b.data[i];
                    
                    //offset2 to left of diagonal
                    offDiag += op1.data[column*(i-off2)+2] * b.data[i-off2];
                    
                    //offset1 to left of diagonal
                    offDiag += op1.data[column*(i-off1)+1] * b.data[i-off1];
                    
                    //offset1 to right of diagonal
                    offDiag += op1.data[column*i+1] * b.data[i+off1];
                    
                    //b(i) = (op2(i) - offDiag(i))/diag(i)
                    //
                    //where offDiag is all the offdiagonal elements of op1 multiplied by corresponding elements of b, diag(i) is leading diagonal element of op1
                    b.data[i] = (op2.data[i] - offDiag)/op1.data[column*i];
                    maxDiff = std::max(maxDiff, double(fabs(b.data[i] - oldb)));
                    
                }
                
                for(int i=size-off1; i<size; i++){
                    
                    offDiag = 0;
                    oldb = b.data[i];
                    
                    //offset2 to left of diagonal
                    offDiag += op1.data[column*(i-off2)+2] * b.data[i-off2];
                    
                    //offset1 to left of diagonal
                    offDiag += op1.data[column*(i-off1)+1] * b.data[i-off1];
                    
                    //b(i) = (op2(i) - offDiag(i))/diag(i)
                    //
                    //where offDiag is all the offdiagonal elements of op1 multiplied by corresponding elements of b, diag(i) is leading diagonal element of op1
                    b.data[i] = (op2.data[i] - offDiag)/op1.data[column*i];
                    
                    maxDiff = std::max(maxDiff, double(fabs(b.data[i] - oldb)));
                    
                }
                
                if (maxDiff < tol){
                    printf("solver exiting after reaching required tolerance %f after %d iterations\n", maxDiff, n);
                    return;
                }
                              
                n++;
                
            }
            break;
            
        case 7:
            
            //loop until error is acceptable or we have exceeded max. iterations
            while (n < maxIter) {
                
                maxDiff = 0;
                for(int i=0; i<off1; i++){
                    
                    offDiag = 0;
                    oldb = b.data[i];
                    
                    //offset1 to right of diagonal
                    offDiag += op1.data[column*i+1] * b.data[i+off1];
                    
                    //offset2 to right of diagonal
                    offDiag += op1.data[column*i+2] * b.data[i+off2];
                    
                    //offset3 to right of diagonal
                    offDiag += op1.data[column*i+3] * b.data[i+off3];
                    
                    //b(i) = (op2(i) - offDiag(i))/diag(i)
                    //
                    //where offDiag is all the offdiagonal elements of op1 multiplied by corresponding elements of b, diag(i) is leading diagonal element of op1
                    b.data[i] = (op2.data[i] - offDiag)/op1.data[column*i];
                    maxDiff = std::max(maxDiff, double(fabs(b.data[i] - oldb)));
                    
                }
                
                for(int i=off1; i<off2; i++){
                    
                    offDiag = 0;
                    oldb = b.data[i];
                    
                    //offset1 to left of diagonal
                    offDiag += op1.data[column*(i-off1)+1] * b.data[i-off1];
                    
                    //offset1 to right of diagonal
                    offDiag += op1.data[column*i+1] * b.data[i+off1];
                    
                    //offset2 to right of diagonal
                    offDiag += op1.data[column*i+2] * b.data[i+off2];
                    
                    //offset3 to right of diagonal
                    offDiag += op1.data[column*i+3] * b.data[i+off3];
                    
                    //b(i) = (op2(i) - offDiag(i))/diag(i)
                    //
                    //where offDiag is all the offdiagonal elements of op1 multiplied by corresponding elements of b, diag(i) is leading diagonal element of op1
                    b.data[i] = (op2.data[i] - offDiag)/op1.data[column*i];
                    maxDiff = std::max(maxDiff, double(fabs(b.data[i] - oldb)));
                    
                }
                
                for(int i=off2; i<off3; i++){
                    
                    offDiag = 0;
                    oldb = b.data[i];
                    
                    //offset2 to left of diagonal
                    offDiag += op1.data[column*(i-off2)+2] * b.data[i-off2];
                    
                    //offset1 to left of diagonal
                    offDiag += op1.data[column*(i-off1)+1] * b.data[i-off1];
                    
                    //offset1 to right of diagonal
                    offDiag += op1.data[column*i+1] * b.data[i+off1];
                    
                    //offset2 to right of diagonal
                    offDiag += op1.data[column*i+2] * b.data[i+off2];
                    
                    //offset3 to right of diagonal
                    offDiag += op1.data[column*i+3] * b.data[i+off3];
                    
                    //b(i) = (op2(i) - offDiag(i))/diag(i)
                    //
                    //where offDiag is all the offdiagonal elements of op1 multiplied by corresponding elements of b, diag(i) is leading diagonal element of op1
                    b.data[i] = (op2.data[i] - offDiag)/op1.data[column*i];
                    maxDiff = std::max(maxDiff, double(fabs(b.data[i] - oldb)));
                    
                }
                
                for(int i=off3; i<size-off3; i++){
                    
                    offDiag = 0;
                    oldb = b.data[i];
                    
                    //offset3 to left of diagonal
                    offDiag += op1.data[column*(i-off3)+3] * b.data[i-off3];
                    
                    //offset2 to left of diagonal
                    offDiag += op1.data[column*(i-off2)+2] * b.data[i-off2];
                    
                    //offset1 to left of diagonal
                    offDiag += op1.data[column*(i-off1)+1] * b.data[i-off1];
                    
                    //offset1 to right of diagonal
                    offDiag += op1.data[column*i+1] * b.data[i+off1];
                    
                    //offset2 to right of diagonal
                    offDiag += op1.data[column*i+2] * b.data[i+off2];
                    
                    //offset3 to right of diagonal
                    offDiag += op1.data[column*i+3] * b.data[i+off3];
                    
                    //b(i) = (op2(i) - offDiag(i))/diag(i)
                    //
                    //where offDiag is all the offdiagonal elements of op1 multiplied by corresponding elements of b, diag(i) is leading diagonal element of op1
                    b.data[i] = (op2.data[i] - offDiag)/op1.data[column*i];
                    maxDiff = std::max(maxDiff, double(fabs(b.data[i] - oldb)));
                    
                }
                
                for(int i=size-off3; i<size-off2; i++){
                    
                    offDiag = 0;
                    oldb = b.data[i];
                    
                    //offset3 to left of diagonal
                    offDiag += op1.data[column*(i-off3)+3] * b.data[i-off3];
                    
                    //offset2 to left of diagonal
                    offDiag += op1.data[column*(i-off2)+2] * b.data[i-off2];
                    
                    //offset1 to left of diagonal
                    offDiag += op1.data[column*(i-off1)+1] * b.data[i-off1];
                    
                    //offset1 to right of diagonal
                    offDiag += op1.data[column*i+1] * b.data[i+off1];
                    
                    //offset2 to right of diagonal
                    offDiag += op1.data[column*i+2] * b.data[i+off2];
                    
                    //b(i) = (op2(i) - offDiag(i))/diag(i)
                    //
                    //where offDiag is all the offdiagonal elements of op1 multiplied by corresponding elements of b, diag(i) is leading diagonal element of op1
                    b.data[i] = (op2.data[i] - offDiag)/op1.data[column*i];
                    maxDiff = std::max(maxDiff, double(fabs(b.data[i] - oldb)));
                    
                }
                
                for(int i=size-off2; i<size-off1; i++){
                    
                    offDiag = 0;
                    oldb = b.data[i];
                    
                    //offset3 to left of diagonal
                    offDiag += op1.data[column*(i-off3)+3] * b.data[i-off3];
                    
                    //offset2 to left of diagonal
                    offDiag += op1.data[column*(i-off2)+2] * b.data[i-off2];
                    
                    //offset1 to left of diagonal
                    offDiag += op1.data[column*(i-off1)+1] * b.data[i-off1];
                    
                    //offset1 to right of diagonal
                    offDiag += op1.data[column*i+1] * b.data[i+off1];
                    
                    //b(i) = (op2(i) - offDiag(i))/diag(i)
                    //
                    //where offDiag is all the offdiagonal elements of op1 multiplied by corresponding elements of b, diag(i) is leading diagonal element of op1
                    b.data[i] = (op2.data[i] - offDiag)/op1.data[column*i];
                    maxDiff = std::max(maxDiff, double(fabs(b.data[i] - oldb)));
                    
                    
                }
                
                for(int i=size-off1; i<size; i++){
                    
                    offDiag = 0;
                    oldb = b.data[i];
                    
                    //offset3 to left of diagonal
                    offDiag += op1.data[column*(i-off3)+3] * b.data[i-off3];
                    
                    //offset2 to left of diagonal
                    offDiag += op1.data[column*(i-off2)+2] * b.data[i-off2];
                    
                    //offset1 to left of diagonal
                    offDiag += op1.data[column*(i-off1)+1] * b.data[i-off1];
                    
                    //b(i) = (op2(i) - offDiag(i))/diag(i)
                    //
                    //where offDiag is all the offdiagonal elements of op1 multiplied by corresponding elements of b, diag(i) is leading diagonal element of op1
                    b.data[i] = (op2.data[i] - offDiag)/op1.data[column*i];
                    maxDiff = std::max(maxDiff, double(fabs(b.data[i] - oldb)));
                    
                }
                
                if (maxDiff < tol){
                    printf("solver exiting after reaching required tolerance %f after %d iterations", maxDiff, n);
                    return;
                }
                
                n++;
                
            }
            break;
            
            
    }
    
    printf("Max number of %d allowed iterations was exceeded, last maxDiff was %f\n", maxIter,maxDiff);
    return;
    
};


template <class T1>
void applyPreconditioner(const symm_band<T1>& A, const diagonal<T1>& Ei, const vector<T1>& r, vector<T1>& z){
  
    int size = A.rows;
    int column = A.columns;
    
    int off1 = A.offset1;
    int off2 = A.offset2;
    int off3 = A.offset3;
    
    vector<T1> q(size);
    
    T1 t;
    
    switch(A.bandwidth){
    
    case 5:
        
        //First solve Lq = r
        //in this range, q(i-1,j) = q(i,j-1) = 0 because they are off-grid
        for(int i=0; i<off1; i++){
            
            if (A.at(i,0) == 0)
                continue;
            
            t = r.data[i];
            
            q.data[i] = t*Ei.data[i];
            
        }
        //in this range, q(i,j-1) = 0 because they are off-grid
        for(int i=off1; i<off2; i++){
            
            if (A.at(i,0) == 0)
                continue;
            
            t = r.data[i] - A.data[(i-off1)*column+1]*Ei.data[i-off1]*q.data[i-off1];
           
            q.data[i] = t*Ei.data[i];
            
        }
        //in this range, all q(i-1,j), q(i,j-1) are on-grid and accessible
        for(int i=off2; i<size; i++){
            
            if (A.at(i,0) == 0)
                continue;
            
            t = r.data[i] - A.data[(i-off1)*column+1]*Ei.data[i-off1]*q.data[i-off1] - A.data[(i-off2)*column+2]*Ei.data[i-off2]*q.data[i-off2];
           
            q.data[i] = t*Ei.data[i];
            
        }
        
        //Next solve L^Tz = q
        for(int i=size-1; i>=size-off1; i--){
            
            if (A.at(i,0) == 0)
                continue;
            
            t = q.data[i];
            
            z.data[i] = t*Ei.data[i];
        
        }
        for(int i=size-off1-1; i>=size-off2; i--){
            
            if (A.at(i,0) == 0)
                continue;
            
            t = q.data[i] - A.data[i*column+1]*Ei.data[i]*z.data[i+off1];
            
            z.data[i] = t*Ei.data[i];
            
        }
        for(int i=size-off2-1; i>=0; i--){
            
            if (A.at(i,0) == 0)
                continue;
            
            t = q.data[i] - A.data[i*column+1]*Ei.data[i]*z.data[i+off1] - A.data[i*column+2]*Ei.data[i]*z.data[i+off2];
            
            z.data[i] = t*Ei.data[i];
            
        }
        break;
    
    case 7:
        
        //First solve Lq = r
        //in this range, q(i-1,j,k) = q(i,j-1,k) = q(i,j,k-1) = 0 because they are off-grid
        for(int i=0; i<off1; i++){
            
            if (A.at(i,0) == 0)
                continue;
            
            t = r.data[i];
            
            q.data[i] = t*Ei.data[i];
            
        }
        //in this range, q(i,j-1,k) = q(i,j,k-1) = 0 because they are off-grid
        for(int i=off1; i<off2; i++){
            
            if (A.at(i,0) == 0)
                continue;
            
            t = r.data[i] - A.data[(i-off1)*column+1]*Ei.data[i-off1]*q.data[i-off1];
            
            q.data[i] = t*Ei.data[i];
            
        }
        //in this range, q(i,j,k-1) = 0 because they are off-grid
        for(int i=off2; i<off3; i++){
            
            if (A.at(i,0) == 0)
                continue;
            
            t = r.data[i] - A.data[(i-off1)*column+1]*Ei.data[i-off1]*q.data[i-off1] - A.data[(i-off2)*column+2]*Ei.data[i-off2]*q.data[i-off2];
            
            q.data[i] = t*Ei.data[i];
            
        }
        //in this range, all q(i,j,k) are on-grid and accessible
        for(int i=off3; i<size; i++){
            
            if (A.at(i,0) == 0)
                continue;
            
            t = r.data[i] - A.data[(i-off1)*column+1]*Ei.data[i-off1]*q.data[i-off1] - A.data[(i-off2)*column+2]*Ei.data[i-off2]*q.data[i-off2]
                          - A.data[(i-off3)*column+3]*Ei.data[i-off3]*q.data[i-off3];
            
            q.data[i] = t*Ei.data[i];
            
        }
        
        
        //Next solve L^Tz = q
        for(int i=size-1; i>=size-off1; i--){
            
            if (A.at(i,0) == 0)
                continue;
            
            t = q.data[i];
            
            z.data[i] = t*Ei.data[i];
            
        }
        for(int i=size-off1-1; i>=size-off2; i--){
            
            if (A.at(i,0) == 0)
                continue;
            
            t = q.data[i] - A.data[i*column+1]*Ei.data[i]*z.data[i+off1];
            
            z.data[i] = t*Ei.data[i];
            
        }
        for(int i=size-off2-1; i>=size-off3; i--){
            
            if (A.at(i,0) == 0)
                continue;
            
            t = q.data[i] - A.data[i*column+1]*Ei.data[i]*z.data[i+off1] - A.data[i*column+2]*Ei.data[i]*z.data[i+off2];
            
            z.data[i] = t*Ei.data[i];
            
        }
        for(int i=size-off3-1; i>=0; i--){
            
            if (A.at(i,0) == 0)
                continue;
            
            t = q.data[i] - A.data[i*column+1]*Ei.data[i]*z.data[i+off1] - A.data[i*column+2]*Ei.data[i]*z.data[i+off2]
                          - A.data[i*column+3]*Ei.data[i]*z.data[i+off3];
            
            z.data[i] = t*Ei.data[i];
            
        }
        break;
            
    };
    
};


template <class T1>
void PCGSolver(const symm_band<T1>& A, const diagonal<T1>& Ei, const vector<T1>& b, vector<T1>& p, T1 tol, int maxIter){
    
    p.fill(0);
    
    vector<T1> r = b;
    
    //r = r - A*p;
    
    vector<T1> z(b.rows);
    applyPreconditioner(A,Ei,r,z);
    
    vector<T1> s = z;
    
    T1 sigma = z*r;
    
    T1 maxr, alpha, beta, sigmanew;
    
    int n = 0;
    
    while (n < maxIter) {
    
        z = A*s;
        
        alpha = sigma/(z*s);
        
        p = p + alpha*s;
        
        r = r - alpha*z;
        
        maxr = r.max();
        
        if (maxr < tol){
            printf("solver exiting after reaching required tolerance %f after %d iterations\n", maxr, n);
            return;
        }
        
        //Apply preconditioner to r, put result in z
        applyPreconditioner(A,Ei,r,z);
        
        sigmanew = z*r;
        
        beta = sigmanew / sigma;
        
        s = z + beta*s;
        
        sigma = sigmanew;
        
        n++;
    
    }
   
    printf("Max number of %d allowed iterations was exceeded, last maxDiff was %f\n", maxIter,maxr);
    return;
    
};


template <class T1>
void CGSolver(const symm_band<T1>& A, const vector<T1>& b, vector<T1>& p, T1 tol, int maxIter){
    
    p.fill(0);
    
    vector<T1> r = b;
    
    //r = r - A*p;
    
    vector<T1> s(b.rows);
    
    vector<T1> z(b.rows);

    s = r;
    
    T1 sigma = r*r;
    
    T1 maxr, alpha, beta, sigmanew;
    
    int n = 0;
    
    while (n < maxIter) {
       
        z = A*s;
        
        alpha = sigma/(s*z);
        
        p = p + alpha*s;
        
        //DEBUG Calculates the residual every frame to minimize cumulative error
        //r = b - A*p;
        
        r = r - alpha*z;
        
        maxr = r.max();
        
        if (maxr < tol){
            printf("solver exiting after reaching required tolerance %f after %d iterations\n", maxr, n);
            return;
        }
        
        sigmanew = r*r;
        
        beta = sigmanew / sigma;
        
        s = r + beta*s;
        
        sigma = sigmanew;
        
        n++;
        
        //DEBUG
        if(n==1){
            
            //Write out r to a binary file on the desktop
            std::ofstream outbin( "/Users/JohnnyK/Desktop/r_binary.bin", std::ios::binary );
            outbin.write( (char *) r.data, sizeof( double )*r.rows );
            outbin.close();
            
            
            //Write out A to a binary file on the desktop
            std::ofstream outbin1( "/Users/JohnnyK/Desktop/A_binary.bin", std::ios::binary );
            outbin1.write( (char *) A.data, sizeof( double )*A.rows*A.columns );
            outbin1.close();
            
        }
        
        
    }
    
    printf("Max number of %d allowed iterations was exceeded, last maxDiff was %f\n", maxIter,maxr);
    return;
    
};


//Scalar multiplication aB = C
template <class T1>
symm_band<T1> operator*(const T1 & scalar, const symm_band<T1>& op2){
 
    symm_band<T1> temp = op2;
    
    for(int i=0; i<temp.rows*temp.columns; i++){
        
        temp.data[i] *= scalar;
        
    }

    return temp;
  
};


template <class T1>
void MIC0precon(const symm_band<T1>& A, diagonal<T1>& Ei){
    
    Ei.fill(0);
    
    double tau = 0.97;
    double sigma = 0.25;
    double e;

    double E1,E2,E3,A11,A22,A33,A12,A13,A21,A23,A31,A32;
    
    int size = A.rows;
    
    int off1 = A.offset1;
    int off2 = A.offset2;
    int off3 = A.offset3;
    
    switch(A.bandwidth){
            
            
        case 5:
            
            //in this range, q(i-1,j) = q(i,j-1) = 0 because they are off-grid
            for(int i=0; i<off1; i++){
                
                e = A.at(i,0);
                
                if (e == 0)
                    continue;
                
                Ei.data[i] = 1.0/sqrt(e);
            }
            
            //in this range, q(i,j-1) = 0 because they are off-grid
            for(int i=off1; i<off2; i++){
                
                e = A.at(i,0);
                
                if (e == 0)
                    continue;
                
                A11 = A.at(i-off1,1);
                A12 = A.at(i-off1,2);
                E1 = Ei.data[i-off1];
                
                e += - A11*A11*E1*E1
                     - tau*(A11*A12*E1*E1);
                
                if (e < sigma*A.at(i,0))
                    e = A.at(i,0);
                Ei.data[i] = 1.0/sqrt(e);
            }
            
            //in this range, all q(i-1,j), q(i,j-1) are on-grid and accessible
            for(int i=off2; i<size; i++){
                
                e = A.at(i,0);
                
                if (e == 0)
                    continue;
                
                A11 = A.at(i-off1,1);
                A12 = A.at(i-off1,2);
                A21 = A.at(i-off2,1);
                A22 = A.at(i-off2,2);
                
                E1 = Ei.data[i-off1];
                E2 = Ei.data[i-off2];
                
                e += - A11*A11*E1*E1
                     - A22*A22*E2*E2
                     - tau*(A11*A12*E1*E1 + A22*A21*E2*E2);
                
                if (e < sigma*A.at(i,0))
                    e = A.at(i,0);
                Ei.data[i] = 1.0/sqrt(e);
            }
            break;
            
    
        case 7:
            
            //in this range, q(i-1,j,k) = q(i,j-1,k) = q(i,j,k-1) = 0 because they are off-grid
            for(int i=0; i<off1; i++){
                
                e = A.at(i,0);
                
                if (e == 0)
                    continue;
                
                Ei.data[i] = 1.0/sqrt(e);
            }
            
            //in this range, q(i,j-1,k) = q(i,j,k-1) = 0 because they are off-grid
            for(int i=off1; i<off2; i++){
                
                e = A.at(i,0);
                
                if (e == 0)
                    continue;
                
                A11 = A.at(i-off1,1);
                A12 = A.at(i-off1,2);
                A13 = A.at(i-off1,3);
                E1 = Ei.data[i-off1];
                
                e += - A11*A11*E1*E1
                     - tau*(A11*(A12 + A13)*E1*E1);
                
                if (e < sigma*A.at(i,0))
                    e = A.at(i,0);
                Ei.data[i] = 1.0/sqrt(e);
            }
            
            //in this range, q(i,j,k-1) = 0 because they are off-grid
            for(int i=off2; i<off3; i++){
                
                e = A.at(i,0);
                
                if (e == 0)
                    continue;
                
                A11 = A.at(i-off1,1);
                A12 = A.at(i-off1,2);
                A13 = A.at(i-off1,3);
                A21 = A.at(i-off2,1);
                A22 = A.at(i-off2,2);
                A23 = A.at(i-off2,3);
                
                E1 = Ei.data[i-off1];
                E2 = Ei.data[i-off2];
                
                e += - A11*A11*E1*E1
                     - A22*A22*E2*E2
                     - tau*(A11*(A12 + A13)*E1*E1 + A22*(A21 + A23)*E2*E2);
                
                if (e < sigma*A.at(i,0))
                    e = A.at(i,0);
                Ei.data[i] = 1.0/sqrt(e);
            }
            
            //in this range, all q(i-1,j,k), q(i,j-1,k) & q(i,j,k-1) are on-grid and accessible
            for(int i=off3; i<size; i++){
                
               
                e = A.at(i,0);
                
                if (e == 0)
                    continue;
                
                A11 = A.at(i-off1,1);
                A12 = A.at(i-off1,2);
                A13 = A.at(i-off1,3);
                A21 = A.at(i-off2,1);
                A22 = A.at(i-off2,2);
                A23 = A.at(i-off2,3);
                A31 = A.at(i-off3,1);
                A32 = A.at(i-off3,2);
                A33 = A.at(i-off3,3);
                
                E1 = Ei.data[i-off1];
                E2 = Ei.data[i-off2];
                E3 = Ei.data[i-off3];
                
                e += - A11*A11*E1*E1
                     - A22*A22*E2*E2
                     - A33*A33*E3*E3
                     - tau*(A11*(A12 + A13)*E1*E1 + A22*(A21 + A23)*E2*E2 + A33*(A31 + A32)*E3*E3);
                
                if (e < sigma*A.at(i,0))
                    e = A.at(i,0);
                Ei.data[i] = 1.0/sqrt(e);
            }
            break;
     
    
    }
    
    
};



// explicit instantiations

template class matrix<int>;
template class matrix<double>;
template matrix<double> operator+<double>(const matrix<double>&, const diagonal<double>&);
template matrix<double> operator*<double>(const matrix<double>&, const diagonal<double>&);
template matrix<double> operator*(const double & scalar, const matrix<double>& op2);
template vector<double> operator*<double>(matrix<double> const&, vector<double> const&);

template class diagonal<int>;
template class diagonal<double>;
template diagonal<double> operator*(const double & scalar, const diagonal<double>& op2);

template void GaussSeidelSolver(const symm_band<double>& op1, const vector<double>& op2, vector<double>& b, double tol, int maxIter);
template void MIC0precon<double>(symm_band<double> const&, diagonal<double>&);
template void PCGSolver<double>(symm_band<double> const&, diagonal<double> const&, vector<double> const&, vector<double>&, double, int);
template void CGSolver<double>(symm_band<double> const&, vector<double> const&, vector<double>&, double, int);

template symm_band<double> operator*(const double & scalar, const symm_band<double>& op2);
template symm_band<double>& symm_band<double>::operator=(const symm_band<double>& op2);
template double symm_band<double>::at(int, int) const;
template double& symm_band<double>::at(int, int);
template symm_band<double>::symm_band(int, int, int, int, int);
template void symm_band<double>::fill(double);
template double symm_band<double>::max();
template bool symm_band<double>::testDiagonal();

template double vector<double>::at(int) const;
template double& vector<double>::at(int);
template vector<double>::vector(int);
template vector<double>::vector(int, double);
template double vector<double>::sum();
template double vector<double>::max();
