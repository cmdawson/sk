#ifndef _INC_SU2_CMP_H
#define _INC_SU2_CMP_H

#include <complex>

#include "matrix.h"
#include "su2.h"

typedef std::complex<double> cdouble;

//! su2_equiv determines whether or not two matrices in SU(2) are equivalent
/*! Two matrices A,B are 'equivalant', up to a prescribed accuracy epsilon
 * if they are within trace distance epsilon of each other */
class su2_equiv 
{
public:
	//! Prescribed accuracy
	double epsilon, gamma;

	//! How close two unitaries have to be before we call them 
	//! Constructor as above
	su2_equiv(double e) : epsilon(e*e), gamma((2.0-e)*(2.0-e)) {}

	//! Comparison operator returns true if and only if A !< B and B!< A
	bool operator()(const Matrix<cdouble>& A, const Matrix<cdouble>& B) const
	{
		double a[] ={real(A(0,0)),-1.0*imag(A(0,1)),real(A(1,0)),imag(A(1,1))};
		double b[] ={real(B(0,0)),-1.0*imag(B(0,1)),real(B(1,0)),imag(B(1,1))};

		a[0] -= b[0];
		a[1] -= b[1];
		a[2] -= b[2];
		a[3] -= b[3];

		double dist = a[0]*a[0] + a[1]*a[1] + a[2]*a[2] + a[3]*a[3];
		if (dist < epsilon || dist > gamma)
			return true;
		else
			return false;

	}
};

#endif
