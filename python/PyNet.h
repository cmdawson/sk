#include <boost/python.hpp>
#include <boost/python/def.hpp>
#include <boost/python/module.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/list.hpp>
#include <boost/python/dict.hpp>

#include <iostream>

#include "Net.h"

using namespace boost::python;
using namespace std;

class PyNet : public Net
{
public:
	typedef complex<double> cdouble;

	//! Usual constructor via a desired tile 'width'
	PyNet(double d) : Net(d) {}

	//! Accepts a filename and loads a pre-saved Net from it
	PyNet(string s) : Net(s.c_str()) {}

	//! DESTROY THEM ALL
	~PyNet() {}

	//! Save a generated Net into a file
	void save(std::string file)
	{
		Net::save(file.c_str());
	}

	//! Accepts a python dictionary mapping gate labels => matrices.
	/*! Matrices should be passed as a flat list */
	void generate(dict& d, int max_depth)
	{
		char msg[128];

		int nm = extract<int>(d.attr("__len__")());

		cdouble** basis = new cdouble*[nm];
		for (int i=0;i<nm;i++)
			basis[i] = new cdouble[4];
		char* labels = new char[nm];

		for (int idx=0;idx<nm;idx++)
		{
			// Get the (label,matrix) tuple
			tuple lm = d.popitem();

			// Store the label
			labels[idx] = extract<char>(lm[0]);

			// Is the 'matrix' passed as a list?
			extract<boost::python::list> ltest(lm[1]);
			if (ltest.check())
			{
				boost::python::list mlist = ltest();

				// check it's got 4 elements
				if (extract<int>(mlist.attr("__len__")()) != 4)
				{
					sprintf(msg,"Object associated with '%c' must be a list"
						" of length 4", labels[idx]); throw Exception(msg);
					return;
				}
				
				basis[idx][0] = extract<cdouble>(mlist[0]);
				basis[idx][1] = extract<cdouble>(mlist[1]);
				basis[idx][2] = extract<cdouble>(mlist[2]);
				basis[idx][3] = extract<cdouble>(mlist[3]);
			}

			else
			{
				sprintf(msg,"Object associated with '%c' must be a list"
					" of length 4", labels[idx]); throw Exception(msg);
				return;
			}
		}

		// Generate the net
		Net::generate(labels,basis,nm,max_depth);

		for (int i=0;i<nm;i++)
			delete[] basis[i];
		delete[] basis;
		delete[] labels;
	}

	//! Return the closest approximating sequence to a matrix 
	/*! The matrix needs to be passed as a flat list of length 4 */
	std::string approximate(boost::python::list& m)
	{
		cdouble elmts[4];

		if (extract<int>(m.attr("__len__")()) != 4)
		{
			char msg[128];
			sprintf(msg,"Net::approximate requires a matrix passed as a "
				" flat list of length 4"); throw Exception(msg);
			return "";
		}

		elmts[0] = extract<cdouble>(m[0]);
		elmts[1] = extract<cdouble>(m[1]);
		elmts[2] = extract<cdouble>(m[2]);
		elmts[3] = extract<cdouble>(m[3]);

		Matrix<cdouble> E(elmts,2,2);
		
		return nearest_string(E);
	}

	//! Evaluate a string and return the resulting matrix as a flat list
	boost::python::list evaluate(std::string s)
	{
		boost::python::list lm;

		Matrix<cdouble> M = Net::evaluate(s);

		lm.append(M(0,0));
		lm.append(M(0,1));
		lm.append(M(1,0));
		lm.append(M(1,1));

		return lm;
	}

	//! Wrap the Solovay-Kitaev function
	/*! Python has no concept of the Net::knot, or of the Matrix<cdouble>, so
	 * we'll return everything as a flat list and deal with it in the sk 
	 * module */

	boost::python::list solovay_kitaev(boost::python::list& m, int d)
	{
		cdouble elmts[4];
		boost::python::list ret;

		if (extract<int>(m.attr("__len__")()) != 4)
		{
			char msg[128];
			sprintf(msg,"Net::approximate requires a matrix passed as a "
				" flat list of length 4"); throw Exception(msg);
			return ret;
		}

		elmts[0] = extract<cdouble>(m[0]);
		elmts[1] = extract<cdouble>(m[1]);
		elmts[2] = extract<cdouble>(m[2]);
		elmts[3] = extract<cdouble>(m[3]);

		Matrix<cdouble> E(elmts,2,2);

		knot ske = Net::solovay_kitaev(E,d);

		ret.append(ske.l);
		ret.append(ske.word);
		ret.append(ske.M(0,0));
		ret.append(ske.M(0,1));
		ret.append(ske.M(1,0));
		ret.append(ske.M(1,1));

		return ret;
	}

};
