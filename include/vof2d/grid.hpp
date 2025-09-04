#ifndef VOF2D_GRID_HPP
#define VOF2D_GRID_HPP

#include <vector>
#include <cstddef>
#include <cmath>
#include <algorithm>

namespace vof2d {

struct GridSize {
    int nx{0}, ny{0};
};

struct Params {
    double Lx{1.0}, Ly{1.0};
    GridSize size{64, 32};
    double dt{1e-3};
    double rho{1000.0};
    double mu{0.0};
    double sigma{0.0};
    double g{9.81};
    bool use_viscosity{false};
    bool use_surface_tension{false};
};

struct Field2D {
    int nx{0}, ny{0};
    std::vector<double> data;

    Field2D() = default;
    Field2D(int nx_, int ny_, double v=0.0): nx(nx_), ny(ny_), data((nx_)*(ny_), v) {}

    inline double& operator()(int i, int j) { return data[j*nx + i]; }
    inline const double& operator()(int i, int j) const { return data[j*nx + i]; }
};

struct MACGrid {
    int nx{0}, ny{0};
    double dx{0.0}, dy{0.0};

    // Staggered velocities
    Field2D u; // size (nx+1, ny)
    Field2D v; // size (nx, ny+1)

    // Cell-centered fields
    Field2D p;   // (nx, ny)
    Field2D c;   // VOF fraction (nx, ny)

    explicit MACGrid(const Params& p_) {
        nx = p_.size.nx;
        ny = p_.size.ny;
        dx = p_.Lx / nx;
        dy = p_.Ly / ny;
        u = Field2D(nx+1, ny, 0.0);
        v = Field2D(nx, ny+1, 0.0);
        p = Field2D(nx, ny, 0.0);
        c = Field2D(nx, ny, 0.0);
    }
};

}

#endif
