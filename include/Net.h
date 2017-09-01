#ifndef _INC_NET_H
#define _INC_NET_H

#include <fstream>
#include <map>
#include <list>
#include <string>
#include <complex>
#include <iostream>

#include "matrix.h"
#include "su2_cmp.h"
#include "su2.h"

#include "Exception.h"

//! Net class
/*! The purpose of a Net is to generate and store instruction sequences from 
 * a given instruction set, and to provide a means for obtaining the nearest 
 * instruction sequence to a given U in SU(2). 
 *
 * A net is initially constructed by specifying a <i>tile width</i>. The 
 * coordinate space for SU(2) is divided up into tiles, so that when we want to
 * look up an approximation to a given U we need only search the tile that it 
 * is in. 
 *
 * When generating a Net, the <i>tile coordinate</i> for each instruction
 * sequence is determined. There are 16 corners to this tile, and the 
 * instruction sequence is associated with each corner. When given a operator
 * U to be approximated we find the closest corner, and search through all the
 * instruction sequences associated with that corner.
 *
 * To generate a Net, we require an instruction set specified as a 2D array
 * containing the instruction gates, and an array of chars containing the 
 * labels we want to denote these matrices by. */
class Net
{
public:
	//! Integer coordinates for an SU(2) matrix in 4D space
	struct icoord
	{
		short x : 16;
		short t : 16;
		short z : 16;
		short y : 16;
	};

	//! Comparison operator for the 4D integer coordinates.
	/*! This uses a kludge to reduce the number of required comparisons. */
	struct icoord_cmp
	{
		bool operator()(const icoord& A, const icoord& B) const
		{
			int* pA = (int *) &A;
			int* pB = (int *) &B;

			if (*pA < *pB)
				return true;
			else if (*pA++ > *pB++)
				return false;
			else if (*pA < *pB)
				return true;
			else 
				return false;
			// If your machine has 64-bit integers comment out the above and
			// replace with this:
			// if(*pA < *pB)
			//		return true;
			// else
			// 		return false;
		}
	};


//public:
	//! A 'knot' is a mildly amusing term for a point in a Net
	struct knot
	{
		int l;				//!< Length of the instruction word
		std::string word;	//!< The instruction word 
		//double coord[4];	//! 4D coordinate of the corresponding matrix
		Matrix<cdouble> M;
	};	

// private:
	//! Storage was chosen somewhat arbitrarily as a list<knot*>
	std::list<knot*> knots;
	typedef std::complex<double> cdouble;

	//! Precision of this net
	double tilewidth;

	//! Number of intervals the coordinate axes are divided into
	int G;

	Matrix<cdouble>* gate_set;	//!< The 'instruction set' of matrices
	int N; 						//!< Size of the instruction set
	char* gate_labels; 			//!< Labels for each of the matrices above
	int* gate_orders;			//!< The order of each gate
	int* gate_inverses;			//!< Index of each gate's inverse

	//! Number of net points searched to arrive at the last approximation
	int search_count;			

	//!< Maps a label to its corresponding index in the above (very oversized)
	char label_indices[512];			

	//! Internal helper for generating nets
	void generate_net(int max_length);

	//! Add the knot (Matrix,word,length) to the net
	void add(Matrix<cdouble>& U, char* word, int l);

//public:

	//! Return the nearest integer coordinate to a Matrix U
	/*! The matrix U lies in a 4D tile bounded by 16 coordinates. If the
	 * optional 'corner' paramemter is given it specifies which of these
	 * 16 coordinates to return. Otherwise it returns the nearest. */
	icoord coordinate(const Matrix<cdouble>& U, int rounding=-1);

	//!< Storage for the net
	typedef std::map<icoord, std::list<knot*> , icoord_cmp> nettype;
	nettype su2net;	

	//! Initialize a Net with a given tile width.
	Net(double tw);

	//! Load a previously generated Net from a file
	Net(const char* file);

	//! DESTROY!
	~Net(void);

	//! Generate the net.
	/*! This function takes an array of n chars which denote the labels
	 * to be used for the net (upper case only). For each label[i] there is
	 * a cdouble vector basis[i] which is the corresponding matrix elements
	 * as u<sub>00</sub>, u<sub>01</sub>, u<sub>10</sub>, u<sub>11</sub>. 
	 * max_length is the maximum length of the instruction sequences we wish
	 * to generate. */
	void generate(char* labels, cdouble** basis, int n, int max_length);

	//! As above, but taking constant size arrays as arguments for convenience
	void generate(char labels[], cdouble basis[][4], int n, int max_length)
	{
		cdouble** b = new cdouble*[n];
		for (int i=0;i<n;i++)
		{
			b[i] = new cdouble[4];
			memcpy(b[i],basis[i],4*sizeof(cdouble));
		}

		generate((char*)labels, b, n, max_length);

		for (int i=0;i<n;i++)
			delete[] b[i];
		delete[] b;
	}

	//! Load a net saved from a file.
	int load(const char* file);

	//! Save the generated net into a file
	void save(const char* file);

	//! Return the nearest knot to U
	knot* nearest(const Matrix<cdouble>& U);

	//! Return the closest approximating string to U in SU(2)
	std::string nearest_string(const Matrix<cdouble>& U)
	{
		return nearest(U)->word;
	}

	//! Evaluate an instruction string and return the resulting Matrix
	Matrix<cdouble> evaluate(const std::string& s);

	//! Invert an instruction string (same as taking its adjoint)
	std::string invert(const std::string& s);

	//! Return the total number of knots seached last approximation
	int searched(void) { return search_count; }

	//! Number of points in the net
	int size(void) { return knots.size(); }

	//! ...
	knot solovay_kitaev(const Matrix<cdouble>& U, int depth);

	//! Output operator
	/*! Prints some terse information about the gate set */
	friend std::ostream& operator<<(std::ostream& out, const Net& N);
};
#endif
