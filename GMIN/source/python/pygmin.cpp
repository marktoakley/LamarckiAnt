//#include "nparray.h"
#include <boost/python.hpp>
#include <math.h>
#include "nparray.h"

using namespace wales;

extern "C" {
    void pygmin_init_();
    double py_get_energy_(double *x); 
    double py_get_energy_gradient_(double *x, double *grad); 
    int py_get_natoms_();
    void pap_lattice_(double *x, double *g, double *energy, int *gtest);
    void py_gmin_get_coords_(double *);
    void userpot_dump_configuration_ (char *str, double * coords, long int strlen);
}

//void potential(boost::python::numeric::array &x, boost::python::numeric::array &grad)

double getEnergy(boost::python::numeric::array &px)
{
    NPArray<1> x(px);
    return py_get_energy_(&x[0]);
}

double getEnergyGradient(boost::python::numeric::array &px, boost::python::numeric::array &pg)
{
    NPArray<1> x(px);
    NPArray<1> g(pg);
    return py_get_energy_gradient_(&x[0], &g[0]);
}

double getPAPEnergy(boost::python::numeric::array &px)
{
    NPArray<1> x(px);
    double energy;
    int l=0;
    pap_lattice_(&x[0], NULL, &energy, &l);
    return energy;
}

int getDOF()
{
  return 3*py_get_natoms_();
}

double getPAPEnergyGradient(boost::python::numeric::array &px, boost::python::numeric::array &pg)
{
    NPArray<1> x(px);
    NPArray<1> g(pg);
    double energy;
    int l=1;
    pap_lattice_(&x[0], &g[0], &energy, &l);
    return energy;
}

void getCoords(boost::python::numeric::array &px)
{
    NPArray<1> x(px);
    py_gmin_get_coords_(&x[0]);
}

void initialize(void)
{
    bool initialized=false;
    if(!initialized)
        pygmin_init_();
    initialized=true;
}

void userpot_dump(const std::string &filename, boost::python::numeric::array &px)
{
    NPArray<1> x(px);
    userpot_dump_configuration_ ((char*)filename.c_str(), &x[0], filename.size());
}

// this is a workaround, is it a old version g++ bug?
#if GMIN_PYTHON_MODULE == 3 
BOOST_PYTHON_MODULE(ambgmin_)
#elif GMIN_PYTHON_MODULE == 2
BOOST_PYTHON_MODULE(dmagmin_)
#elif GMIN_PYTHON_MODULE == 4 
BOOST_PYTHON_MODULE(oxdnagmin_)
#else
BOOST_PYTHON_MODULE(gmin_)
#endif
{
  using namespace boost::python;
  import_array();
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
  def("initialize", initialize);
  def("getEnergy", getEnergy);
  def("getEnergyGradient", getEnergyGradient);
  def("getNAtoms", py_get_natoms_);
  def("getDOF", getDOF);
  def("getPAPEnergy", getPAPEnergy);
  def("getPAPEnergyGradient", getPAPEnergyGradient);
  def("getCoords", getCoords);
  def("userpot_dump", userpot_dump);
}
