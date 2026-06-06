#include "grid.h"


bool CGSolver(const SymmBandMatrix& A, const SymmBandMatrix::VectorType& b, SymmBandMatrix::VectorType& p, double tol, int maxIter);
const SymmBandMatrix::VectorType operator+(const SymmBandMatrix::VectorType& lhs, const SymmBandMatrix::VectorType& rhs);
const SymmBandMatrix::VectorType operator-(const SymmBandMatrix::VectorType& lhs, const SymmBandMatrix::VectorType& rhs);
const SymmBandMatrix::VectorType operator*(const double& lhs, const SymmBandMatrix::VectorType& rhs);
const double max(const SymmBandMatrix::VectorType& x);