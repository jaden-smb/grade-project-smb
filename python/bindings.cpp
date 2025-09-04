#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include "vof2d/solver.hpp"

namespace py = pybind11;
using namespace vof2d;

PYBIND11_MODULE(vof2d, m) {
    m.doc() = "2D VOF incompressible Navier–Stokes solver (educational)";

    py::class_<GridSize>(m, "GridSize")
        .def(py::init<int,int>(), py::arg("nx"), py::arg("ny"))
        .def_readwrite("nx", &GridSize::nx)
        .def_readwrite("ny", &GridSize::ny);

    py::class_<Params>(m, "Params")
        .def(py::init<>())
        .def_readwrite("Lx", &Params::Lx)
        .def_readwrite("Ly", &Params::Ly)
        .def_readwrite("size", &Params::size)
        .def_readwrite("dt", &Params::dt)
        .def_readwrite("rho", &Params::rho)
        .def_readwrite("mu", &Params::mu)
        .def_readwrite("sigma", &Params::sigma)
        .def_readwrite("g", &Params::g)
        .def_readwrite("use_viscosity", &Params::use_viscosity)
        .def_readwrite("use_surface_tension", &Params::use_surface_tension);

    py::class_<Field2D>(m, "Field2D")
        .def_property_readonly("nx", [](const Field2D& f){ return f.nx; })
        .def_property_readonly("ny", [](const Field2D& f){ return f.ny; })
        .def("to_numpy", [](const Field2D& f){
            return py::array_t<double>({f.ny, f.nx}, {sizeof(double)*f.nx, sizeof(double)}, f.data.data());
        }, py::return_value_policy::reference_internal);

    py::class_<MACGrid>(m, "MACGrid")
        .def_property_readonly("nx", [](const MACGrid& g){ return g.nx; })
        .def_property_readonly("ny", [](const MACGrid& g){ return g.ny; })
        .def_property_readonly("dx", [](const MACGrid& g){ return g.dx; })
        .def_property_readonly("dy", [](const MACGrid& g){ return g.dy; })
        .def_property_readonly("u", [](MACGrid& g) -> Field2D& { return g.u; }, py::return_value_policy::reference_internal)
        .def_property_readonly("v", [](MACGrid& g) -> Field2D& { return g.v; }, py::return_value_policy::reference_internal)
        .def_property_readonly("p", [](MACGrid& g) -> Field2D& { return g.p; }, py::return_value_policy::reference_internal)
        .def_property_readonly("c", [](MACGrid& g) -> Field2D& { return g.c; }, py::return_value_policy::reference_internal);

    py::class_<Solver2D>(m, "Solver2D")
        .def(py::init<const Params&>())
        .def("step", &Solver2D::step)
        .def("compute_mass", &Solver2D::compute_mass)
        .def_property_readonly("grid", [](Solver2D& s) -> MACGrid& { return s.grid(); }, py::return_value_policy::reference_internal)
        .def("initialize_vof", [](Solver2D& s, py::function f){
            s.initialize_vof([&](double x, double y){ return f(x,y).cast<double>(); });
        })
        .def("initialize_velocity", [](Solver2D& s, py::function fu, py::function fv){
            s.initialize_velocity([&](double x, double y){ return fu(x,y).cast<double>(); },
                                  [&](double x, double y){ return fv(x,y).cast<double>(); });
        });
}
