#include <iostream>

#include "Net.h"
#include "Exception.h"

using namespace std;
typedef complex<double> cdouble;

int main(void)
{
	// Square root of -1
	cdouble I = cdouble(0,1);

	// Set up the labels array
	char labels[] = {'H', 'T'};

	// And set up the 2D array of matrices correspnding to these labels
	cdouble gates[][4] =
	{
	  { I/SQRT2, I/SQRT2, I/SQRT2, -1.0*I/SQRT2 }, 	
	  { cdouble(cos(PI/8),-sin(PI/8)) ,0, 0, cdouble(cos(PI/8),sin(PI/8)) } 
	};

	// Construct a Net with the HT gate set and a tile width of 0.14
	Net enet(0.18);
	cout << "... generating net ..." << endl;
	enet.generate(labels,gates,2,14);
	
	// Optionally save the net into a file
	// enet.save("ht14.net"); 
	// Which can be loaded next time as follows
	// Net enet("lps7.net");

	cout << enet << endl;

	// Generate a random 4D unit vector
	srand(time(0));
	double a[4];
	double norm = 0;
	for (int i=0;i<4;i++) {a[i] = (float)rand()/RAND_MAX; norm += a[i]*a[i]; }
	norm = sqrtf(norm);
	for (int i=0;i<4;i++) { a[i] /= norm; }

	// and construct a matrix from this unit vector
	Matrix<cdouble> A(2,2);
	su2::cart4_to_mat(a,A);

	cout << "A = " << A << endl;

	try {
	
	cout << "calling Solovay-Kitaev() to depth 5" << endl;
	Net::knot ska = enet.solovay_kitaev(A,5);

	cout << "got a sequence of length " << ska.l << endl;

	cout << "A ~ " << ska.M << endl;
	cout << "accuracy " << su2::proj_trace_dist(A,ska.M) << endl;

	// Confirm this sequence is what it's supposed to be
	// Matrix<cdouble> B = enet.evaluate(ska.word);
	// cout << "0 = " << su2::proj_trace_dist(ska.M,B) << endl;
	
	// catch any exceptions
	}
	catch (Exception& E) { cout << "*** ERROR: " << E.msg << endl; }
}
