//#include "nparray.h"
#include <boost/python.hpp>
#include <math.h>
#include "nparray.h"
#include <string>

using namespace wales;

extern "C" {
    void pygmin_init_();
    double py_get_energy_(double *x); 
    double py_get_energy_gradient_(double *x, double *grad); 
    int py_get_natoms_();
    int py_get_dof_();
    void py_dmagmin_get_coords_(double *);
    void py_dmagmin_takestep_(double *x);
    void dmacrys_reduce_cell_rigid_(double *x, int *was_changed);

    void __genrigid_MOD_transformrigidtoc(int *min, int *max, double *x, double *xrigid);
    void __genrigid_MOD_transformctorigid(double *x, double *xrigid);

    extern int __genrigid_MOD_nrigidbody;
    void __dmacrys_interface_MOD_dmacrys_dump_cif(char *filename, double * coords, char *cifname, long int strlen1, long int strlen2);
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

void getCoords(boost::python::numeric::array &px)
{
    NPArray<1> x(px);
    py_dmagmin_get_coords_(&x[0]);
}

void takestep(boost::python::numeric::array &px)
{
    NPArray<1> x(px);
    py_dmagmin_takestep_(&x[0]);
}

void initialize(void)
{
    static bool initialized=false;
    if(!initialized) {
        pygmin_init_();
	//py_dmamgin_init();
    }
    
    initialized=true;
}

bool reduceCell(boost::python::numeric::array &px)
{
    NPArray<1> x(px);
    int changed;
    
    dmacrys_reduce_cell_rigid_(&x[0], &changed);
    return changed!=0;
}

int getNRigidBody() {
    return __genrigid_MOD_nrigidbody;
}

void toRigid(boost::python::numeric::array &px, boost::python::numeric::array &pxrigid)
{
    NPArray<1> x(px);
    NPArray<1> xrigid(pxrigid);
    
    __genrigid_MOD_transformctorigid(&x[0], &xrigid[0]);
}

void toAtomistic(boost::python::numeric::array &px, boost::python::numeric::array &pxrigid)
{
    NPArray<1> x(px);
    NPArray<1> xrigid(pxrigid);
    
    int cmin=1, cmax=getNRigidBody();
    __genrigid_MOD_transformrigidtoc(&cmin, &cmax, &x[0], &xrigid[0]);
}

void writeCIF(const std::string &filename, boost::python::numeric::array &px, const std::string &cifname)
{
    NPArray<1> x(px);
    __dmacrys_interface_MOD_dmacrys_dump_cif((char*)filename.c_str(), &x[0], (char*)cifname.c_str(), filename.size(), cifname.size());
}

BOOST_PYTHON_MODULE(dmagmin_)
{
  using namespace boost::python;
  import_array();
  boost::python::numeric::array::set_module_and_type("numpy", "ndarray");
  def("initialize", initialize);
  def("getEnergy", getEnergy);
  def("getEnergyGradient", getEnergyGradient);
  def("getNAtoms", py_get_natoms_);
  def("getDOF", py_get_dof_);
  def("getCoords", getCoords);
  def("takestep", takestep);
  def("getNRigidBody", getNRigidBody);
  def("reduceCell", reduceCell);
  def("toAtomistic", toAtomistic);
  def("toRigid", toRigid);
  def("writeCIF", writeCIF);
}
