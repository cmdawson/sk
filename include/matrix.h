/* A Quantum Information Matrix Toolkit
 * 
 * The last thing the world needs is another C++ Matrix library. However ...
 *
 * This code is intended to facilitate coding C++ numerics related to Quantum
 * Information. Naturally you'd probably be better off using Matlab or   
 * something but I like C++ and enjoy coding in it. This library provides:
 * 	- A templated Matrix class with all the usual arithmetic operators and
 * some additional functions commonly used in Quantum Mechanics (such as the
 * partial trace and tensor product etc).
 *	- A convenient interface to certain Lapack linear algebra routines for
 * inverses, eigen and singular value decompositions.
 *	- Routines for the calculation of various entanglement measures. 
 *
 * Source code:
 *    http://www.physics.uq.edu.au/people/dawson/matrix/matrix-0.2.tar.gz
 *
 * Makefiles are provided for Linux (tested with g++ 2.95 and 3.2.1) and Windows
 * (Visual C++ 7) which compile shared objects / dlls. You'll also need the 
 * Lapack and Blas libraries which are provided by most major linux
 * distributions. For Windows I recommend the CLaw32 stuff from
 * http://www.netlib.org/clapack/
 *
 * Alternatively for windows you can download precompiled dlls
 *
 * Win32 libraries:
 *    http://www.physics.uq.edu.au/people/dawson/matrix/matrixw32-0.2.zip
 * 
 * This document is somewhat incomplete.  
 *
 * If you make any improvements, fix any bugs, or anything please let me know.
 *
 * dawson@physics.uq.edu.au */

#ifndef _MATRIX_
#define _MATRIX_
#include <complex>
#include <iostream>
#include <fstream>
#include <cstdlib>

#define PI 3.1415926535897932384626
#define SQRT2 1.41421356237309504880169
#define EABS(a) (abs(a)<1e-15?0.0:abs(a))
#define EVAL(a) (abs(a)<1e-14?0.0:a)
//#define SWAP(a,b) {temp=a;a=b;b=temp;}

template <typename Z> static Z gconj(const Z &other) { return other; }
template <> static std::complex<double>
gconj<std::complex<double> >(const std::complex<double> & other) {
	 return std::conj(other);
}

/*! Matrix Class 
 * Takes any built in type or standard library types such as
 * std::complex<double> as its template argument.  */

template <typename Z> class Matrix {
public:

	Z *A;				//!< Pointer to the matrix elements 
	unsigned int M;		//!< Number of Rows
	unsigned int N;		//!< Number of Columns

	Z* scratch;			//!< Scratch space for working

	Matrix(void) : M(0), N(0), A(0), scratch(0) {}

	//! Empty Constructor
	/*! Takes two arguments specifying the rows and columns */ 
	Matrix(unsigned int Rows, unsigned int Cols);

	//! Initilizing constructor
	/*! This takes an additional pointer to an initilizing array containing
	 * the matrix elements. For example \code int a[] = {1,2,3,4}; 
	 * Matrix<int> A(a,2,2); \endcode */
	Matrix(Z *elements ,unsigned int Rows, unsigned int Cols);

	Matrix(const Matrix &A);				//!< Assignment Constructor
	Matrix& operator =(const Matrix &C);	//!< Copy Constructor

	~Matrix();

	//! Elementwise Access
	inline Z&  operator() (unsigned int i, unsigned int j);
	const Z&  operator() (unsigned i, unsigned j) const { return A[j*M+i]; }

	//! In place addition. 
	/*! This performs A = A+B. It is usually preferably to do addition like
 	 * this when you're worried about performance as it doesn't require any
  	 * temporary variables or copying of return values */ 
	Matrix& operator+=(const Matrix& M);

	//! In place subtraction (see above)
	Matrix& operator-=(const Matrix& M);

	//! In place multiplication
	/* Again this is the more efficient way to do multiplication. This 
	 * should take a const reference. Unfortunately this clashes with the
	 * overloads below and I have no idea how to reconcidle them. */
	Matrix& operator*=(Matrix M);

	//! In place scalar multiplication.
	/* This is templated so you can conveniently multiply a std::complex matrix
	 * by an integer and so on. (Naturally you wouldn't want to try it the
	 * other way round) */

	#ifdef _MSC_VER
	template<typename S>
	Matrix& operator*=(const S &c) {
		for (unsigned i=0;i<M*N;i++)
			A[i] *= c;
		return *this;
	}
	#else
	template <typename S> Matrix& operator*=(const S& s);
	#endif

	//! Output operator
	#ifdef _MSC_VER
	friend std::ostream& operator<<(std::ostream&, const Matrix<Z>&);
	#else
	friend std::ostream& operator<<<>(std::ostream&, const Matrix&);
	#endif

	//! Dump matrix to filename	
	void dump(const char* filename);	
	//! Read matrix from file created by dump(const char*)
	void read(const char* filename);
};

/*! \file matrix.h
 * \brief Basic Matrix functionality.
 * Provides arithmetic operators, various matrix functions and an output 
 * operator. */

/*! @name Arithmetic Operators
 * These are all the basic matrix arithmetic operations that you'd expect.
 * Note that the templated scalar multiplication makes actual matrix 
 * multiplication a partial specialization. This might be too much for older
 * MS compilers. */
//@{
template <typename Z>
Matrix<Z> operator+(const Matrix<Z>&, const Matrix<Z>&);

template <typename Z>
Matrix<Z> operator-(const Matrix<Z> &, const Matrix<Z> &);


template <typename Z>
Matrix<Z> operator*(const Z &, const Matrix<Z> &);

#ifndef _MSC_VER
template <typename Z, typename S>
Matrix<Z> operator*(const S&, const Matrix<Z> &);
#else
template <typename Z, typename S>
Matrix<Z> operator*(const S &a, const Matrix<Z> &B) {
	Matrix<Z> R = B;
	return R*=a;
}
#endif

template <typename Z>
Matrix<Z> operator*(const Matrix<Z> &, const Matrix<Z> &);
//@}

/*! @name Matrix Functions
 * Basic matrix functions like the Adjoint, Determinant etc that don't require
 * any involved calculation. For the 'other' sort of matrix function see the
 * Lapack interface
 * Note that all these functions have reference arguments. This is to avoid
 * an expensive copying of the Matrix argument, but contains a little gotcha.
 * On most compilers (MS VC++ being the exception) this means that things like
 * \code t = Trace(Adjoint(M)) \endcode aren't possible. You should do 
 * \code AM = Adjoint(M); t = Trace(AM); \endcode instead. This wouldn't happen
 * if the functions took const refererences, but then we wouldn be able to use
 * the element access operator M(i,j) in these functions. 'Unresolved Issue' :P
 */
//@{

//! Adjoint or Transpose (depending on Z) 
/*! For std::complex matrices this calculates the adjoint, for real or integer
 * matrices it returns the Transpose */
template <typename Z>
Matrix<Z> Adjoint(const Matrix<Z>&);

//! Normal Trace 
template <typename Z>
Z Trace(Matrix<Z>&);

/*! \brief Partial Trace
 *
 * This is a dumb but sound way of computing the partial trace 
 * over B of an ABC system from
 *   \f[\sum_i (I_A \otimes <i| \otimes I_C) R (I_A \otimes |i> \otimes I_C)\f]
 * dA represents the dimension of the first subsystem. dB the 
 * dimension of the system being traced over. If there is any left
 * over, this is considered subsystem C */
template <typename Z>
Matrix<Z> PTrace(Matrix<Z>& R, unsigned dA, unsigned dB);

//! Tensor (or Kronecker or Outer) Product \f$A \otimes B\f$
template <typename Z>
Matrix<Z> TensorProduct(Matrix<Z>& A, Matrix<Z>& B);

/*! \brief Determinant
 *
 * Recursively calculates the determinant by expanding along a row. This
 * scales horribly with matrix dimension. */ 
template <typename Z>
Z Determinant(Matrix<Z>& M); 

/*! \brief Submatrix formed by removing a given row and column 
 *
 * The Matrix S is formed from M by removing row i and column j */
template <typename Z>
void SubMatrix(Matrix<Z>& M, Matrix<Z>& S, unsigned i, unsigned j);
//@}

/*! @name Output Operator */
//@{
//! Output to any ostream (normally cout).  
//template <typename Z>
//std::ostream& operator<<(std::ostream& output, const Matrix<Z>& A);
//@}

/* ----------------------------------------------------------------------- *
 * End of Prototypes and DOXYGEN comments
 * Start of Implementation
 * ----------------------------------------------------------------------- */

template<typename Z>
Matrix<Z>::Matrix(unsigned int m, unsigned int n) : M(m), N(n) {
	A = new Z[m*n];
	scratch = new Z[n];
}

template<typename Z>
Matrix<Z>::Matrix(Z *a,unsigned int m, unsigned int n) : M(m), N(n) {
	A = new Z[m*n];
	scratch = new Z[n];
	for (unsigned i=0;i<n;i++) 
		for (unsigned j=0;j<m;j++)
			A[i*m+j] = a[n*j+i];
}

template<typename Z>
Matrix<Z>::Matrix(const Matrix<Z> &t) : M(t.M), N(t.N) {
	A = new Z[t.M*t.N];
	scratch = new Z[N];
	for (unsigned i=0;i<t.M*t.N;i++) A[i] = t.A[i];
}

template<typename Z>
Matrix<Z>& Matrix<Z>::operator =(const Matrix<Z> &t) {
	if (this != &t) {
		if (N != t.N || M != t.M) {
			delete[] A;
			delete[] scratch;
			A = new Z[t.M*t.N];
			scratch = new Z[t.N];
			M = t.M; N = t.N;
		}
		for (unsigned i=0;i<M*N;i++) A[i] = t.A[i];
	}
	return *this;
}

template<typename Z>
Matrix<Z>::~Matrix() {
	delete[] A;
	delete[] scratch;
} 

template<typename Z>
inline Z& Matrix<Z>::operator() (unsigned int i, unsigned int j) {
	return A[j*M+i];
}

template<typename Z>
Matrix<Z>& Matrix<Z>::operator+=(const Matrix<Z> &B) {
	for (unsigned i=0;i<M*N;i++)
		A[i] += B.A[i];
	return *this;
}

template<typename Z>
Matrix<Z>& Matrix<Z>::operator-=(const Matrix<Z> &B) {
	for (unsigned i=0;i<M*N;i++)
		A[i] -= B.A[i];
	return *this;
}

template<typename Z>
Matrix<Z>& Matrix<Z>::operator*=(Matrix<Z> B) {

	// First check this makes sense
	if (N != B.M)
		std::cerr << "Incompatible matrices passed to Matrix *=\n";

	// And ensure we're not post-multiplying by ourself
	if (this == &B)
	{
		Matrix<Z> TMP = B;
		return operator*=(TMP);
	}

	// If B has more columns than A we'll need to resize
	if (B.N > N)
	{
		Z *tmp = new Z[M*B.N];
		memcpy(tmp,A,M*N*sizeof(Z));
		delete[] A;
		A = tmp;

		delete[] scratch;
		scratch = new Z[B.N];
	}

	for (unsigned i=0;i<M;i++)
	{
		// copy this row of *this  into scratch
		for (unsigned j=0;j<N;j++)
			scratch[j] = A[j*M+i];

		// now for each jth column of B
		for (unsigned j=0;j<B.N;j++) 
		{
			// Set A(i,j) equal to A(i,*) . B(*,j)
			A[i+M*j] = 0;
			for (unsigned k=0;k<N;k++) 
				A[i+M*j] += (scratch[k] * B.A[k+j*B.M]);
		}
	}

	N = B.N;
	return *this;
}

#ifndef _MSC_VER
template<typename Z> template<typename S>
Matrix<Z>& Matrix<Z>::operator*=(const S &c)
{
	for (unsigned i=0;i<M*N;i++)
		A[i] *= c;
	return *this;
}
#endif

template<typename Z> void Matrix<Z>::dump(const char *fn) {
	std::ofstream out_file;

	out_file.open(fn, std::ios::out | std::ios::binary);

	if (out_file.bad()) {
		std::cerr << "Unable to open matrix file for writing" << fn << std::endl;
		exit(0);
	}

	out_file.write(static_cast<unsigned int *>(&M),sizeof(unsigned int));
	out_file.write(static_cast<unsigned int *>(&N),sizeof(unsigned int));

	out_file.write(A,N*M*sizeof(Z));

	if (out_file.bad()) {
		std::cerr << "Unable to write matrix file " << fn << std::endl;
		exit(0);
	}

	out_file.close();
}


template<typename Z> void Matrix<Z>::read(const char *fn) {
	unsigned int m,n;
	std::ifstream in_file;

	in_file.open(fn, std::ios::in | std::ios::binary);

	if (in_file.bad()) {
		std::cerr << "Unable to open matrix file " << fn  
			<< " for reading\n";
		exit(0);
	}

	in_file.read(static_cast<unsigned int *>(&m),sizeof(unsigned int));
	in_file.read(static_cast<unsigned int *>(&n),sizeof(unsigned int));

	// resize if necessary	
	if (n<N || m<M) {
		delete[] A;
		A = new Z[n*m];
	}
	N = n; M = m;

	in_file.read(A,N*M*sizeof(Z));

	if (in_file.bad()) {
		std::cerr << "Unable to read matrix file " << fn << std::endl;
		exit(0);
	}

	in_file.close();
}

/* ----------------------------------------------------------------------- *
 * End of Member Functions
 * ----------------------------------------------------------------------- */

template <typename Z>
Matrix<Z> operator*(const Matrix<Z> &A, const Matrix<Z> &B)
{
	Matrix<Z> R = A;
	return R*=B;
}

template <typename Z>
Matrix<Z> operator+(const Matrix<Z> &A, const Matrix<Z> &B) {
	Matrix<Z> R = A;
	return R+=B;
}

template <typename Z>
Matrix<Z> operator-(const Matrix<Z> &A, const Matrix<Z> &B) {
	Matrix<Z> R = A;
	return R-=B;
}

template <typename Z>
Matrix<Z> operator*(const Z &a, const Matrix<Z> &B) {
	Matrix<Z> R = B;
	//std::cout << "*Z &a\n";
	return R*=a;
}

#ifndef _MSC_VER
template <typename Z, typename S>
Matrix<Z> operator*(const S &a, const Matrix<Z> &B) {
	Matrix<Z> R = B;
	return R*=a;
}
#endif

template <typename Z>
Matrix<Z> Adjoint(const Matrix<Z>& U) {
	Matrix<Z> ADJ(U.N,U.M);
	for (unsigned i=0;i<U.M;i++) 
	    for(unsigned j=0;j<U.N;j++) 
		ADJ(j,i) = gconj<Z>(U(i,j));
	return ADJ;
}

template<typename Z>
Z Trace(Matrix<Z>& A) {
	Z result = 0;
	for (unsigned i=0;i<A.M;i++)
		result += A(i,i);

	return result;
}

template <typename Z>
std::ostream& operator<<(std::ostream& output, const Matrix<Z>& A) {
	output << "[\n";
	for (unsigned i=0;i<A.M;i++) {
		for (unsigned j=0;j<A.N;j++) 
			//output << "\t" << A(i,j) << "\t";
			output << "\t" << A.A[j*A.M+i] << "\t";
		output << "\n";
	}
	output << "]";

	return output; 
}

template <typename Z>
void SubMatrix(Matrix<Z> &P, Matrix<Z> &C, unsigned r, unsigned c) {

	for (unsigned i=0;i<r;i++) {
		for (unsigned j=0;j<c;j++)
			C(i,j) = P(i,j);
		for (unsigned j=c+1;j<P.N;j++)
			C(i,j-1) = P(i,j);
	}

	for (unsigned i=r+1;i<P.M;i++) {
		for (unsigned j=0;j<c;j++)
			C(i-1,j) = P(i,j);
		for (unsigned j=c+1;j<P.N;j++)
			C(i-1,j-1) = P(i,j);
	}
}

template<typename Z>
Matrix<Z> TensorProduct(Matrix<Z>& A,Matrix<Z>& B) {

	if (A.M == 1 && A.N == 1) {
		return B;
	}

	if (B.M == 1 && B.N == 1)
		return A;

	unsigned rows = A.M*B.M;
	unsigned cols = A.N*B.N;

	Matrix<Z> M(rows,cols);

	for (unsigned i=0;i<rows;i++) 
		for (unsigned j=0;j<cols;j++) 
			M(i,j) = A(i/B.M,j/B.N) * B(i%B.M,j%B.N);

	return M;
}


template<typename Z>
Matrix<Z> PTrace(Matrix<Z>& M, unsigned dA, unsigned dB) {
	// check this is actually possible
	//assert((M.M % (dA*dB)) == 0);

	// figure out how big the returned matrix will be
	unsigned D = M.M / dB;
	unsigned dC;

	// if there is an offset of 0, just consider this a
	// 1d system for simplicity
	if (dA == 0)
		dA = 1;

	// figure out the size of subsystem C
	dC = D / dA;

	Matrix<Z> R(D,D), IA(dA,dA), IC(dC,dC), B(1,dB);
	// Tensor product I_A x <e| I_C
	Matrix<Z> T(D,M.M);

	B(0,0) = 1;
	for (unsigned i=1;i<dB;i++)
		B(0,i) = 0;

	for (unsigned i=0;i<dA;i++) {
		IA(i,i) = 1;
		for (unsigned j=0;j<i;j++)
			IA(i,j) = IA(j,i) = 0;
	}

	for (unsigned i=0;i<dC;i++) {
		IC(i,i) = 1;
		for (unsigned j=0;j<i;j++)
			IC(i,j) = IC(j,i) = 0;
	}
		
	// now do the partial trace
	Matrix<Z> TMP = TensorProduct(B,IC);
	T = TensorProduct(IA,TMP);
	R = T*M*Adjoint(T);

	for (unsigned i=1;i<dB;i++) {
		B(0,i-1) = 0; B(0,i) = 1;
		TMP = TensorProduct(B,IC);
		T = TensorProduct(IA,TMP);
		R += (T*M*Adjoint(T));
	}

	return R;
}

template <typename Z>
Z Determinant(Matrix<Z>& M) {
	Z d; int s;
	if (M.M == 2)
		return M(0,0)*M(1,1)-M(1,0)*M(0,1);

	d = (Z)0; s = 1;

	Matrix<Z> S(M.M-1,M.N-1);

	for (unsigned i=0;i<M.M;i++) {
		SubMatrix(M,S,0,i);
		d += static_cast<const Z&>(s)*M(0,i) * Determinant(S);
		s *= -1;
		
	}
	return d;
}

#endif
