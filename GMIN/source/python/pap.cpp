#include <boost/python.hpp>
#include <math.h>
#include "nparray.h"

using namespace wales;

extern "C" {
    void pap_lattice_(double *x, double *g, double *energy, int gtest);
}

double getEnergy(boost::python::numeric::array &px)
{
    NPArray<1> x(px);
    return pap_lattice(&x[0], 0);
}

double getEnergyGradient(boost::python::numeric::array &px, boost::python::numeric::array &pg)
{
    NPArray<1> x(px);
    NPArray<1> g(pg);
    return py_get_energy_gradient_(&x[0], &g[0]);
}

BOOST_PYTHON_MODULE(pygmin)
{
  using namespace boost::python;
  import_array();
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
  def("init", pygmin_init_);
  def("getEnergy", getEnergy);
  def("getEnergyGradient", getEnergyGradient);
  def("getNAtoms", py_get_natoms_);
}

