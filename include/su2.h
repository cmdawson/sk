#ifndef _INC_SU2_H
#define _INC_SU2_H

#include "matrix.h"
#include "complex"

typedef std::complex<double> cdouble;

//! A collection of SU(2) utilities
/*! This provides:
 * 	- a bunch of routines for converting SU(2) matrices to and from various
 * 	parameterizations
 * 	- a bunch of decompositions including the important 'group commutator'
 * 	decomposition.
 *
 * 	The parameterizations are the 3D Bloch vector <b>r</b> from 
 * 	\f$U = e^{i \mathbf{r} . \mathbf{\sigma}}\f$  and the 4D unit vector 
 * 	from $U = aI + ibX + icY + idZ$. The naming convention is pretty self
 * 	explanatory e.g. mat_to_cart3(Matrix<cdouble>& U, double* r) will
 * 	calculate the cartesian form of the Block vector for U and place it in r,
 * 	while polar4_to_mat(double* p, Matrix<cdouble>& U) will construct the
 * 	matrix U from the polar form of a 4D unit vector. */
namespace su2 
{

//! 4D polar coordinates. Not really used.
struct polar { double r, t1, t2, t3; };

//! Calculate cartesian form of the 3D Bloch vector for U
void mat_to_cart3(const Matrix<cdouble>& U, double* c);

//! Construct a matrix U from the 3D cartesian Bloch vector
void cart3_to_mat(double* c, Matrix<cdouble>& U);

//! Calculate polar form of the 3D Bloch vector for U
void mat_to_polar3(const Matrix<cdouble>& U, double* p); 

//! Construct a matrix U from the 3D polar Bloch vector
void polar3_to_mat(double* p, Matrix<cdouble>& U);

//! Convert a polar 3D vector to a cartesian 3D vector
void polar_to_cart3(double* p, double* c);

//! Convert a cartesian 3D vector to a polar 3D vector
void cart_to_polar3(double* c, double* p);

//! Given a 3D cartesian vector z, find vectors x,y s.t. \f$z = x \times y\f$
void x_factor(double *z, double *x, double *y);

//! Calculate cartesian form of the 4D unit vector for U
void mat_to_cart4(const Matrix<cdouble>& U, double* a);

//! Construct a matrix U from the 4D cartesian unit vector
void cart4_to_mat(double *a, Matrix<cdouble>& U);

//! Calculate the polar form of the 4D unit vector for U
void mat_to_polar4(const Matrix<cdouble>&, double *); 

//! Construct a matrix U from the 4D polar unit vector
void polar4_to_mat(double*, Matrix<cdouble>&);

//! Convert a polar form 4D vector to cartesian form
void polar_to_cart4(double*, double*);

//! Convert a cartesian form 4D vector to polar form
void cart_to_polar4(double*, double*);

//! Legacy conversion from a 4D polar struct to cartesian form
void polar_to_cart4(polar& p, double* r);

//! Legacy conversion from a 4D cartesian vector to a polar struct
void polar4_to_mat(polar&, Matrix<cdouble>&);

//! Quick and dirty calculation of the eigenvalues of an SU(2) matrix
void eig_vals(Matrix<cdouble>& U, cdouble *evals);

//! Quick and dirty calculation of the eigenvectors of an SU(2) matrix
void eig_vecs(Matrix<cdouble>& U, cdouble *evecs);

//! Trace distance between SU(2) matrices U,V
double trace_dist(const Matrix<cdouble>& U, const Matrix<cdouble>& V);

//! Trace distance between SU(2) matrices U,V
/*! Except V in this case is given by its 4D cartesian coordinates */
double trace_dist(const Matrix<cdouble>& U, double* coord4);

//! 'Projective' trace distance between SU(2) matrices U,V
double proj_trace_dist(const Matrix<cdouble>& U, const Matrix<cdouble>& V);

//! 'Projective' trace distance between SU(2) matrices U,V
/*! V in this case is given by its 4D cartesian coordinates */
double proj_trace_dist(const Matrix<cdouble>& U, double* coord4);

//! Find 'balanced' V,W in SU(2) such that U = [V,W]
void group_factor(const Matrix<cdouble>& U, Matrix<cdouble>& V, Matrix<cdouble>& W);

//! 'Balanced' group commutator decomposition for U a rotation about the x-axis
void x_group_factor(const Matrix<cdouble>&, Matrix<cdouble>&, Matrix<cdouble>&);

//! Find an SU(2) matrix S such that \f$ V = S W S^\dagger \f$
/*! (or it might be the other way around) */
void similarity_matrix(const Matrix<cdouble>& V, Matrix<cdouble>& W, Matrix<cdouble>&);

//! Norm of a real vector v of length n
double norm(double* v, int n);

};

#endif
