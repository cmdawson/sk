#include <ctime>
#include <cassert>

#include "su2.h"

#ifndef PI
#define PI 3.1415926535897932384626
#endif

#ifndef SQRT2
#define SQRT2 1.41421356237309504880169
#endif

using namespace su2;
using namespace std;

cdouble _basis[] = {1,0,0,1,0,1,1,0,0,cdouble(0,-1),cdouble(0,1),0,1,0,0,-1};

Matrix<cdouble> Id2 = Matrix<cdouble>(_basis,2,2);
Matrix<cdouble> PX = Matrix<cdouble>(_basis+4,2,2);
Matrix<cdouble> PY = Matrix<cdouble>(_basis+8,2,2);
Matrix<cdouble> PZ = Matrix<cdouble>(_basis+12,2,2);

// Find balanced V,W such that U = [V,W]
void su2::group_factor(const Matrix<cdouble>& U, Matrix<cdouble>& V, Matrix<cdouble>& W)
{
	// U is a rotation about some axis, find out by how much, then make
	// that a rotation about the X-axis.
	Matrix<cdouble> XU(2,2), S(2,2), aS(2,2), A(2,2), B(2,2);
	double u[3], n;
	mat_to_cart3(U,u);
	n = norm(u,3);

	double a[] = {n,0,0};
	cart3_to_mat(a,XU);
	similarity_matrix(U,XU,S);
	aS = Adjoint(S);
	x_group_factor(XU,A,B);
	
	V = S*A*aS;
	W = S*B*aS;
}

// This is an experiment using the method of decomposition in Ky Fan  
/*void su2::GroupFactor3(Matrix<cdouble>& U, Matrix<cdouble>& V, Matrix<cdouble>& W)
{
	Matrix<cdouble> aE(2,2), aV(2,2), Y(2,2), C(2,2);
	cdouble v[4], l[2];
	double x[3];
	mat_to_cart3(U,x);
	cout << x[0] << " " << x[1] << " " << x[2] << endl;

	eigVecs(U,v);
	Matrix<cdouble> E(v,2,2);
	aE = Adjoint(E);

	eig_vals(U,l);
	cdouble a[] = {l[0],0,0,l[0]*l[1]};
	cdouble b[] = {conj(l[0]*l[1]),0,0,conj(l[0])};
	Matrix<cdouble> A(a,2,2);
	Matrix<cdouble> B(b,2,2);

	*V = E*A*aE;
	mat_to_cart3(*V,x);
	cout << "A is \n" << x[0] << " " << x[1] << " " << x[2] << endl;
	aV = Adjoint(V);
	Y = E*B*aE;
	mat_to_cart3(aV,x);
	cout << "A' is \n" << x[0] << " " << x[1] << " " << x[2] << endl;

	mat_to_cart3(Y,x);
	cout << "B is " << x[0] << " " << x[1] << " " << x[2] << endl;

	similarity_matrix(Y,aV,W);
}*/

void su2::similarity_matrix(const Matrix<cdouble>& A,Matrix<cdouble>& B,Matrix<cdouble>& S)
{
	double a[3], b[3], s[3];
	double na,nb,ns,ab;
	mat_to_cart3(A,a);
	mat_to_cart3(B,b);

	na = norm(a,3);
	nb = norm(b,3);
	assert (abs(na- nb) < 1e-12);

	ab = 0;
	for (int i=0;i<3;i++) 
		ab += a[i]*b[i];

	s[0] = b[1]*a[2] - a[1]*b[2];
	s[1] = a[0]*b[2] - b[0]*a[2];
	s[2] = b[0]*a[1] - a[0]*b[1];	

	ns = norm(s,3);
	if (abs(ns) < 1e-12) {
		// The vectors are parallel or anti parallel 
		// i am just pretending they are parallel so fix this.
		S = Id2;
		return;
	}
		
	for (int i=0;i<3;i++) 
		s[i] *= (acos(ab/(na*nb))/ns);

	cart3_to_mat(s,S);
} 


void su2::x_group_factor(const Matrix<cdouble>& A, Matrix<cdouble>& B, Matrix<cdouble>& C)
{
	double a[4], b[3], c[3], st, ct, theta, alpha, phi;
	mat_to_cart4(A,a);
	Matrix<cdouble> W(2,2),aB(2,2); 
	
	// st = pow(0.5 - 0.5*sqrt(1 - a[1]*a[1]),0.25);
	st = pow(0.5 - 0.5*a[0],0.25);
	ct = sqrt(1-st*st);
	theta = 2*asin(st);
	alpha = atan(st);

	b[0] = c[0] = theta*st*cos(alpha);
	b[1] = c[1] = theta*st*sin(alpha);	
	b[2] = theta*ct;
	c[2] = -b[2];

	cart3_to_mat(c,B);
	cart3_to_mat(b,W);

	aB = Adjoint(B);

	similarity_matrix(W,aB,C);
}

// Functions to convert matrices into 3 vector representations exp(-ix.sigma/2)
// and vice versa. Also 3-vector polar and cartesian conversions
void su2::mat_to_cart3(const Matrix<cdouble>& U,double *x)
{
	double p, sx1, sx2, sx3, th, sinth, costh;

	sx1 = -1 * imag(U(0,1));
	sx2 = real(U(1,0));
	sx3 = imag((U(1,1) - U(0,0))/2.0);

	costh = real((U(0,0) + U(1,1))/2.0);

	sinth = sqrt(sx1*sx1 + sx2*sx2 + sx3*sx3);

	if (sinth == 0)
	{
		x[0] = 2*acos(costh);
		x[1] = x[2] = 0;
	}
	else
	{
		th = atan2(sinth,costh);
		x[0] = 2*th*sx1/sinth;
		x[1] = 2*th*sx2/sinth;
		x[2] = 2*th*sx3/sinth;
	}
}

void su2::cart3_to_mat(double *x, Matrix<cdouble>& U)
{
	double th = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);

	if (th == 0) 
		U = Id2;
	else
	{
		U = cos(th/2)*Id2 + (-cdouble(0,1)*sin(th/2)*x[0]/th)*PX
				+ (-cdouble(0,1)*sin(th/2)*x[1]/th)*PY
				+ (-cdouble(0,1)*sin(th/2)*x[2]/th)*PZ;
	}

}

void su2::polar3_to_mat(double *p, Matrix<cdouble>& U)
{
	double c[3];
	polar_to_cart3(p,c);
	cart3_to_mat(c,U);
}
	
void su2::mat_to_polar3(const Matrix<cdouble>& U, double *p)
{
	double c[3];
	mat_to_cart3(U,c);
	cart_to_polar3(c,p);
}

void su2::polar_to_cart3(double *polar, double *cart)
{
	double x = polar[0] * sin(polar[2]);
	cart[0] = x * cos(polar[1]);
	cart[1] = x * sin(polar[1]);
	cart[2] = polar[0] * cos(polar[2]);
}

void su2::cart_to_polar3(double *cart, double *polar)
{
	polar[0] = sqrt(cart[0]*cart[0]+cart[1]*cart[1]+cart[2]*cart[2]);
	polar[1] = atan2(cart[1], cart[0]);
	polar[2] = atan2(sqrt(cart[0]*cart[0] + cart[1]*cart[1]), cart[2]);
	
	if (polar[2]<0)
	{
		polar[2] = polar[2] * -1;	
		polar[1] = (polar[1]<PI) ? polar[1]+PI : polar[1]-PI;
	}

	polar[1] = (polar[1]<0) ? polar[1]+(2*PI) : polar[1];	
}


// likewise matrix to 4 vector routines and 4-vector polar/cart stuff
void su2::polar_to_cart4(double *polar, double *cart)
{
	cart[0] = polar[0] * sin(polar[1])*sin(polar[2])*sin(polar[3]);
	cart[1] = polar[0] * cos(polar[1])*sin(polar[2])*sin(polar[3]);
	cart[2] = polar[0] * cos(polar[2])*sin(polar[3]);
	cart[3] = polar[0] * cos(polar[3]);
}

void su2::polar_to_cart4(polar& p, double *cart)
{
	cart[0] = p.r * sin(p.t1)*sin(p.t2)*sin(p.t3);
	cart[1] = p.r * cos(p.t1)*sin(p.t2)*sin(p.t3);
	cart[2] = p.r * cos(p.t2)*sin(p.t3);
	cart[3] = p.r * cos(p.t3);
}

void su2::cart4_to_mat(double *c, Matrix<cdouble>& U)
{
	U = c[0]*Id2 + (c[1]*cdouble(0,-1))*PX + (c[2]*cdouble(0,-1))*PY
		 + (c[3]*cdouble(0,-1))*PZ;
}

void su2::mat_to_cart4(const Matrix<cdouble>& U, double *q)
{
	q[0] = real(U(0,0));
	q[1] = -1 * imag(U(0,1));
	q[2] = real(U(1,0));
	q[3] = imag(U(1,1));

	// dodgy rounding off to 0
	for (int i=0;i<4;i++)
		q[i] = abs(q[i]) < 1e-15 ? 0 : q[i];
}

void su2::cart_to_polar4(double *cart, double *polar) {
	double foo;
	polar[0] = norm(cart,4);

	foo = (polar[0]>0) ? cart[3]/polar[0] : 0;
	foo = (foo<-1) ? -1 : (foo>1) ? 1 : foo;
	polar[3] = acos(foo);

	foo = (polar[0]>0 && sin(polar[3])>0) ? cart[2]/(sin(polar[3])*polar[0]) :0;
	foo = (foo<-1) ? -1 : (foo>1) ? 1 : foo;
	polar[2] = acos(foo);

	polar[1] = atan2(cart[0],cart[1]);
	polar[1] = (polar[1] < 0) ? polar[1]+2*PI : polar[1];
}

void su2::mat_to_polar4(const Matrix<cdouble>& M, double *p)
{
	double c[4];
	mat_to_cart4(M,c);
	cart_to_polar4(c,p);
}

void su2::polar4_to_mat(double *polar, Matrix<cdouble>& M)
{
	double cart[4];
	polar_to_cart4(polar,cart);
	cart4_to_mat(cart,M);
}

void su2::polar4_to_mat(polar& p, Matrix<cdouble>& M)
{
	double c[4];
	polar_to_cart4(p,c);
	cart4_to_mat(c,M);
}

/* ----------------------------------------------------------------------------
 * Miscellaneous stuff */

double su2::norm(double *v, int d)
{
	double n=0;
	for (int i=0;i<d;i++)
		n += v[i]*v[i];

	return sqrt(n);
}
		
void su2::x_factor(double *z, double *x, double *y)
{
	double zp[3],xp[3],yp[3];

	cart_to_polar3(z,zp);
	yp[0] = xp[0] = sqrt(zp[0]);
	yp[1] = PI/2;
	yp[2] = zp[2] + PI/2;
	xp[1] = zp[1] + PI/2;
	xp[2] = zp[2];

	polar_to_cart3(xp, x);
	polar_to_cart3(yp, y);
}

// These next two functions will need to be replaced by some LU  decomposition
// to work for SU(N)

void su2::eig_vals(Matrix<cdouble>& A, cdouble *l)
{
	l[0] = 0.5 * (A(0,0)+A(1,1)+sqrt(pow(A(0,0)-A(1,1),2)+4.0*A(0,1)*A(1,0)));
	l[1] = 0.5 * (A(0,0)+A(1,1)-sqrt(pow(A(0,0)-A(1,1),2)+4.0*A(0,1)*A(1,0)));
}

void su2::eig_vecs(Matrix<cdouble>& A, cdouble *v)
{
	cdouble n1,n2;
	cdouble l[2];
	eig_vals(A,l);
	v[0] = v[1] = 1;
	v[2] = (l[0] - A(0,0))/A(0,1);
	v[3] = (l[1] - A(0,0))/A(0,1);

	n1 = sqrt(1.0 + v[2]*conj(v[2]));
	v[0] /= n1;
	v[2] /= n1;

	n2 = sqrt(1.0 + v[3]*conj(v[3]));
	v[1] /= n2;
	v[3] /= n2;
}

double su2::trace_dist(const Matrix<cdouble>& A, const Matrix<cdouble>& B)
{
	cdouble tr[2];
	Matrix<cdouble> C(2,2), D(2,2);
	C = A-B;
	D = Adjoint(C);

	C*=D;

	eig_vals(C,tr);
	//std::cout << imag(tr[0]) << " " << imag(tr[1]) << "\n";
	assert(abs(imag(tr[0])) < 1e-15 && abs(imag(tr[1])) < 1e-15);
		
	return 0.5 * sqrt(real(*tr)) + sqrt(real(*(tr+1)));
}

double su2::trace_dist(const Matrix<cdouble>& A, double* coord4)
{
	double d = 0.0;

	d += (real(A(0,0))-coord4[0])*(real(A(0,0))-coord4[0]);
	d += (imag(A(0,1))+coord4[1])*(imag(A(0,1))+coord4[1]);
	d += (real(A(1,0))-coord4[2])*(real(A(1,0))-coord4[2]);
	d += (imag(A(1,1))-coord4[3])*(imag(A(1,1))-coord4[3]);

	return sqrtf(d);
}

double su2::proj_trace_dist(const Matrix<cdouble>& A, const Matrix<cdouble>& B)
{
	double d = 0.0;
	double n = 0.0;

	d += (real(A(0,0))-real(B(0,0)))*(real(A(0,0))-real(B(0,0)));
	d += (imag(A(0,1))-imag(B(0,1)))*(imag(A(0,1))-imag(B(0,1)));
	d += (real(A(1,0))-real(B(1,0)))*(real(A(1,0))-real(B(1,0)));
	d += (imag(A(1,1))-imag(B(1,1)))*(imag(A(1,1))-imag(B(1,1)));
			
	n += (real(A(0,0))+real(B(0,0)))*(real(A(0,0))+real(B(0,0)));
	n += (imag(A(0,1))+imag(B(0,1)))*(imag(A(0,1))+imag(B(0,1)));
	n += (real(A(1,0))+real(B(1,0)))*(real(A(1,0))+real(B(1,0)));
	n += (imag(A(1,1))+imag(B(1,1)))*(imag(A(1,1))+imag(B(1,1)));

	return (d < n) ? sqrtf(d) : sqrtf(n);
}


double su2::proj_trace_dist(const Matrix<cdouble>& A, double* coord4)
{
	double d = 0.0;
	double n = 0.0;

	d += (real(A(0,0))-coord4[0])*(real(A(0,0))-coord4[0]);
	d += (imag(A(0,1))+coord4[1])*(imag(A(0,1))+coord4[1]);
	d += (real(A(1,0))-coord4[2])*(real(A(1,0))-coord4[2]);
	d += (imag(A(1,1))-coord4[3])*(imag(A(1,1))-coord4[3]);

	n += (real(A(0,0))+coord4[0])*(real(A(0,0))+coord4[0]);
	n += (imag(A(0,1))-coord4[1])*(imag(A(0,1))-coord4[1]);
	n += (real(A(1,0))+coord4[2])*(real(A(1,0))+coord4[2]);
	n += (imag(A(1,1))+coord4[3])*(imag(A(1,1))+coord4[3]);

	return (d < n) ? sqrtf(d) : sqrtf(n);
}
