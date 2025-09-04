#ifndef VOF2D_SOLVER_HPP
#define VOF2D_SOLVER_HPP

#include "grid.hpp"
#include <vector>

namespace vof2d {

class Solver2D {
public:
    explicit Solver2D(const Params& params);

    const Params& params() const { return params_; }
    MACGrid& grid() { return grid_; }
    const MACGrid& grid() const { return grid_; }

    // Main step: advection + body forces + viscosity (optional) + projection + VOF update
    void step();

    // Initialize VOF fraction using a user callback f(x,y) in [0,1]
    template <typename F>
    void initialize_vof(F&& frac_fun);

    // Set velocities from callback
    template <typename FU, typename FV>
    void initialize_velocity(FU&& ufun, FV&& vfun);

    // Diagnostics
    double compute_mass() const; // sum(c) * dx*dy

private:
    Params params_;
    MACGrid grid_;

    // internals
    void apply_body_forces();
    void apply_viscosity();
    void advect_velocities();
    void project_incompressible();
    void advect_vof_conservative();
    void plic_reconstruct_and_correct();
    void apply_surface_tension();
};

// Template method definitions inline
template <typename F>
inline void Solver2D::initialize_vof(F&& frac_fun) {
    for (int j = 0; j < grid_.ny; ++j) {
        for (int i = 0; i < grid_.nx; ++i) {
            double x = (i+0.5)*grid_.dx;
            double y = (j+0.5)*grid_.dy;
            grid_.c(i,j) = std::clamp(frac_fun(x,y), 0.0, 1.0);
        }
    }
}

template <typename FU, typename FV>
inline void Solver2D::initialize_velocity(FU&& ufun, FV&& vfun) {
    for (int j = 0; j < grid_.ny; ++j) {
        for (int i = 0; i < grid_.nx+1; ++i) {
            double x = i*grid_.dx;
            double y = (j+0.5)*grid_.dy;
            grid_.u(i,j) = ufun(x,y);
        }
    }
    for (int j = 0; j < grid_.ny+1; ++j) {
        for (int i = 0; i < grid_.nx; ++i) {
            double x = (i+0.5)*grid_.dx;
            double y = j*grid_.dy;
            grid_.v(i,j) = vfun(x,y);
        }
    }
}

}

#endif
