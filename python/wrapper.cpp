#include "PyNet.h"
#include "su2.h"

#include <boost/python/exception_translator.hpp>

using namespace boost::python;
using namespace std;

void excpt_translate(const Exception& E)
{
	PyErr_SetString(PyExc_UserWarning, E.msg.c_str());
}


BOOST_PYTHON_MODULE(sk_base)
{
	register_exception_translator<Exception>(&excpt_translate);

	class_<PyNet>("net", init<double>())
		.def(init<std::string>())
		.def("save",&PyNet::save)
		.def("generate_base",&PyNet::generate)
		.def("approximate_base",&PyNet::approximate)
		.def("evaluate_base",&PyNet::evaluate)
		.def("invert",&PyNet::invert)
		.def("searched",&PyNet::searched)
		.def("size",&PyNet::size)
		.def("solovay_kitaev_base",&PyNet::solovay_kitaev)
		// This next one was messing up last I tried it
		//.def(boost::python::self_ns::str(self))
	;

}


