#include "vof2d/solver.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace vof2d {

Solver2D::Solver2D(const Params& params)
    : params_(params), grid_(params) {}

void Solver2D::apply_body_forces() {
    // Gravity in -y on v-component
    const double dt = params_.dt;
    const double g = params_.g;
    auto& v = grid_.v;
    for (int j = 0; j < grid_.ny+1; ++j) {
        for (int i = 0; i < grid_.nx; ++i) {
            v(i,j) += -g * dt; // simplistic body force
        }
    }
}

void Solver2D::apply_viscosity() {
    if (!params_.use_viscosity || params_.mu <= 0.0) return;
    // Simple explicit Laplacian on u and v
    const double dt = params_.dt;
    const double mu = params_.mu;
    const double dx2 = grid_.dx * grid_.dx;
    const double dy2 = grid_.dy * grid_.dy;

    Field2D un = grid_.u;
    Field2D vn = grid_.v;

    // u: (nx+1, ny)
    for (int j = 1; j < grid_.ny-1; ++j) {
        for (int i = 1; i < grid_.nx; ++i) {
            double lap = (un(i+1,j) - 2*un(i,j) + un(i-1,j)) / dx2
                       + (un(i,j+1) - 2*un(i,j) + un(i,j-1)) / dy2;
            grid_.u(i,j) += dt * mu * lap;
        }
    }
    // v: (nx, ny+1)
    for (int j = 1; j < grid_.ny; ++j) {
        for (int i = 1; i < grid_.nx-1; ++i) {
            double lap = (vn(i+1,j) - 2*vn(i,j) + vn(i-1,j)) / dx2
                       + (vn(i,j+1) - 2*vn(i,j) + vn(i,j-1)) / dy2;
            grid_.v(i,j) += dt * mu * lap;
        }
    }
}

void Solver2D::advect_velocities() {
    // Semi-Lagrangian (stable but not strictly conservative). Simpler for demo.
    // For boundedness, clamp interpolations.
    const double dt = params_.dt;
    Field2D u0 = grid_.u;
    Field2D v0 = grid_.v;

    auto sample_u = [&](double x, double y){
        // x on u-grid: i in [0..nx], y in [0..ny-1]
        double gx = x / grid_.dx;
        double gy = y / grid_.dy;
        int i = std::clamp((int)std::floor(gx), 0, grid_.nx);
        int j = std::clamp((int)std::floor(gy), 0, grid_.ny-1);
        int i1 = std::min(i+1, grid_.nx);
        int j1 = std::min(j+1, grid_.ny-1);
        double tx = std::clamp(gx - i, 0.0, 1.0);
        double ty = std::clamp(gy - j, 0.0, 1.0);
        double v00 = u0(i,j), v10 = u0(i1,j), v01 = u0(i,j1), v11 = u0(i1,j1);
        return (1-tx)*((1-ty)*v00 + ty*v01) + tx*((1-ty)*v10 + ty*v11);
    };
    auto sample_v = [&](double x, double y){
        double gx = x / grid_.dx;
        double gy = y / grid_.dy;
        int i = std::clamp((int)std::floor(gx), 0, grid_.nx-1);
        int j = std::clamp((int)std::floor(gy), 0, grid_.ny);
        int i1 = std::min(i+1, grid_.nx-1);
        int j1 = std::min(j+1, grid_.ny);
        double tx = std::clamp(gx - i, 0.0, 1.0);
        double ty = std::clamp(gy - j, 0.0, 1.0);
        double v00 = v0(i,j), v10 = v0(i1,j), v01 = v0(i,j1), v11 = v0(i1,j1);
        return (1-tx)*((1-ty)*v00 + ty*v01) + tx*((1-ty)*v10 + ty*v11);
    };

    // Advect u
    for (int j = 0; j < grid_.ny; ++j) {
        for (int i = 0; i < grid_.nx+1; ++i) {
            double x = i * grid_.dx; // u at face centers: (i*dx, (j+0.5)*dy)
            double y = (j + 0.5) * grid_.dy;
            double ux = sample_u(x, y);
            double vy = sample_v(x, y);
            double xb = x - dt * ux;
            double yb = y - dt * vy;
            grid_.u(i,j) = sample_u(std::clamp(xb, 0.0, grid_.nx * grid_.dx),
                                     std::clamp(yb, 0.0, grid_.ny * grid_.dy));
        }
    }
    // Advect v
    for (int j = 0; j < grid_.ny+1; ++j) {
        for (int i = 0; i < grid_.nx; ++i) {
            double x = (i + 0.5) * grid_.dx;
            double y = j * grid_.dy;
            double ux = sample_u(x, y);
            double vy = sample_v(x, y);
            double xb = x - dt * ux;
            double yb = y - dt * vy;
            grid_.v(i,j) = sample_v(std::clamp(xb, 0.0, grid_.nx * grid_.dx),
                                     std::clamp(yb, 0.0, grid_.ny * grid_.dy));
        }
    }
}

void Solver2D::project_incompressible() {
    // Simple Jacobi pressure solve to enforce divergence-free
    const int iters = 80;
    const double dx = grid_.dx, dy = grid_.dy;
    const double dx2 = dx*dx, dy2 = dy*dy;
    const double invBeta = 1.0 / (2.0*(dx2+dy2));

    // compute divergence
    Field2D div(grid_.nx, grid_.ny, 0.0);
    for (int j = 0; j < grid_.ny; ++j) {
        for (int i = 0; i < grid_.nx; ++i) {
            double du = grid_.u(i+1,j) - grid_.u(i,j);
            double dv = grid_.v(i,j+1) - grid_.v(i,j);
            div(i,j) = (du/dx + dv/dy);
        }
    }

    Field2D p = grid_.p; // initial guess
    for (int k = 0; k < iters; ++k) {
        Field2D pn = p;
        for (int j = 1; j < grid_.ny-1; ++j) {
            for (int i = 1; i < grid_.nx-1; ++i) {
                p(i,j) = ((pn(i+1,j)+pn(i-1,j))*dy2 + (pn(i,j+1)+pn(i,j-1))*dx2 - div(i,j)*dx2*dy2) * invBeta;
            }
        }
        // Neumann at boundaries (dp/dn=0)
        for (int i = 0; i < grid_.nx; ++i) { p(i,0) = p(i,1); p(i,grid_.ny-1) = p(i,grid_.ny-2); }
        for (int j = 0; j < grid_.ny; ++j) { p(0,j) = p(1,j); p(grid_.nx-1,j) = p(grid_.nx-2,j); }
    }
    grid_.p = p;

    // Subtract pressure gradient
    for (int j = 0; j < grid_.ny; ++j) {
        for (int i = 1; i < grid_.nx; ++i) {
            grid_.u(i,j) -= (p(i,j) - p(i-1,j)) / dx;
        }
    }
    for (int j = 1; j < grid_.ny; ++j) {
        for (int i = 0; i < grid_.nx; ++i) {
            grid_.v(i,j) -= (p(i,j) - p(i,j-1)) / dy;
        }
    }
    // Velocity boundary conditions: no-flow across domain boundary
    for (int j = 0; j < grid_.ny; ++j) { grid_.u(0,j) = 0.0; grid_.u(grid_.nx,j) = 0.0; }
    for (int i = 0; i < grid_.nx; ++i) { grid_.v(i,0) = 0.0; grid_.v(i,grid_.ny) = 0.0; }
}

void Solver2D::advect_vof_conservative() {
    // Donor-acceptor with conservative flux limiting and bounding.
    const double dt = params_.dt;
    const double dx = grid_.dx, dy = grid_.dy;

    Field2D cn = grid_.c;

    // Raw face fluid fluxes (signed): amount of fluid volume crossing face in dt
    Field2D Fx(grid_.nx+1, grid_.ny, 0.0);
    for (int j = 0; j < grid_.ny; ++j) {
        for (int i = 0; i < grid_.nx+1; ++i) {
            double uf = grid_.u(i,j);
            int iup = (uf > 0) ? std::max(i-1,0) : std::min(i, grid_.nx-1);
            double c_up = cn(iup,j);
            Fx(i,j) = uf * dt * dy * c_up;
        }
    }
    Field2D Fy(grid_.nx, grid_.ny+1, 0.0);
    for (int j = 0; j < grid_.ny+1; ++j) {
        for (int i = 0; i < grid_.nx; ++i) {
            double vf = grid_.v(i,j);
            int jup = (vf > 0) ? std::max(j-1,0) : std::min(j, grid_.ny-1);
            double c_up = cn(i,jup);
            Fy(i,j) = vf * dt * dx * c_up;
        }
    }

    // Face scaling factors from outflow/inflow constraints
    Field2D Sx(grid_.nx+1, grid_.ny, 1.0);
    Field2D Sy(grid_.nx, grid_.ny+1, 1.0);
    const double vol = dx*dy;

    for (int j = 0; j < grid_.ny; ++j) {
        for (int i = 0; i < grid_.nx; ++i) {
            double cvol = cn(i,j) * vol;
            // outgoing contributions
            double outL = std::max(-Fx(i,j), 0.0);
            double outR = std::max( Fx(i+1,j), 0.0);
            double outB = std::max(-Fy(i,j), 0.0);
            double outT = std::max( Fy(i,j+1), 0.0);
            double sumOut = outL + outR + outB + outT;
            if (sumOut > cvol && sumOut > 0.0) {
                double s = cvol / sumOut;
                if (outL > 0) Sx(i,j) = std::min(Sx(i,j), s);
                if (outR > 0) Sx(i+1,j) = std::min(Sx(i+1,j), s);
                if (outB > 0) Sy(i,j) = std::min(Sy(i,j), s);
                if (outT > 0) Sy(i,j+1) = std::min(Sy(i,j+1), s);
            }

            // incoming capacity limit
            double vcap = (1.0 - cn(i,j)) * vol;
            double inL = std::max( Fx(i,j), 0.0);
            double inR = std::max(-Fx(i+1,j), 0.0);
            double inB = std::max( Fy(i,j), 0.0);
            double inT = std::max(-Fy(i,j+1), 0.0);
            double sumIn = inL + inR + inB + inT;
            if (sumIn > vcap && sumIn > 0.0) {
                double s = vcap / sumIn;
                if (inL > 0) Sx(i,j) = std::min(Sx(i,j), s);
                if (inR > 0) Sx(i+1,j) = std::min(Sx(i+1,j), s);
                if (inB > 0) Sy(i,j) = std::min(Sy(i,j), s);
                if (inT > 0) Sy(i,j+1) = std::min(Sy(i,j+1), s);
            }
        }
    }

    // Apply scaling
    for (int j = 0; j < grid_.ny; ++j) {
        for (int i = 0; i < grid_.nx+1; ++i) Fx(i,j) *= Sx(i,j);
    }
    for (int j = 0; j < grid_.ny+1; ++j) {
        for (int i = 0; i < grid_.nx; ++i) Fy(i,j) *= Sy(i,j);
    }

    // Update cells
    for (int j = 0; j < grid_.ny; ++j) {
        for (int i = 0; i < grid_.nx; ++i) {
            double inflow = Fx(i,j) - Fx(i+1,j) + Fy(i,j) - Fy(i,j+1);
            double cnew = cn(i,j) + inflow / vol;
            grid_.c(i,j) = std::clamp(cnew, 0.0, 1.0);
        }
    }
}

void Solver2D::plic_reconstruct_and_correct() {
    // Placeholder: in a full PLIC, reconstruct interface plane per cell using gradients
    // and correct fluxes to preserve boundedness and sharpness. For now, we perform
    // a light sharpening limited diffusion when 0<c<1 to counteract smearing.
    Field2D cn = grid_.c;
    Field2D c2 = cn;
    const double alpha = 0.01; // small sharpening factor
    for (int j = 1; j < grid_.ny-1; ++j) {
        for (int i = 1; i < grid_.nx-1; ++i) {
            double lap = (cn(i+1,j) + cn(i-1,j) + cn(i,j+1) + cn(i,j-1) - 4*cn(i,j));
            double cc = cn(i,j) + alpha * lap;
            c2(i,j) = std::clamp(cc, 0.0, 1.0);
        }
    }
    grid_.c = c2;
}

void Solver2D::apply_surface_tension() {
    if (!params_.use_surface_tension || params_.sigma <= 0.0) return;
    // Continuum Surface Force (CSF): F = sigma * kappa * grad(c)
    // Estimate normals n = grad(c)/|grad(c)| and curvature kappa = div(n)
    const int nx = grid_.nx, ny = grid_.ny;
    Field2D nxF(nx, ny, 0.0), nyF(nx, ny, 0.0);
    Field2D kappa(nx, ny, 0.0);
    const double dx = grid_.dx, dy = grid_.dy;
    // gradients
    for (int j = 1; j < ny-1; ++j) {
        for (int i = 1; i < nx-1; ++i) {
            double dcx = (grid_.c(i+1,j) - grid_.c(i-1,j)) / (2*dx);
            double dcy = (grid_.c(i,j+1) - grid_.c(i,j-1)) / (2*dy);
            double nrm = std::sqrt(dcx*dcx + dcy*dcy) + 1e-12;
            nxF(i,j) = dcx / nrm;
            nyF(i,j) = dcy / nrm;
        }
    }
    // curvature div(n)
    for (int j = 1; j < ny-1; ++j) {
        for (int i = 1; i < nx-1; ++i) {
            double dnxdx = (nxF(i+1,j) - nxF(i-1,j)) / (2*dx);
            double dnydy = (nyF(i,j+1) - nyF(i,j-1)) / (2*dy);
            kappa(i,j) = dnxdx + dnydy;
        }
    }
    // Add force to velocities (explicit): u += dt * (F_x / rho), v += dt * (F_y / rho)
    const double dt = params_.dt;
    const double rho = params_.rho;
    const double sigma = params_.sigma;
    // Interpolate F at faces using central differences of c scaled by curvature
    // Fx ~ sigma * kappa * dcdx, Fy ~ sigma * kappa * dcdy
    for (int j = 0; j < ny; ++j) {
        for (int i = 1; i < nx; ++i) {
            int ic = std::min(std::max(i,1), nx-1);
            int jc = std::min(std::max(j,1), ny-1);
            double dcdx = (grid_.c(ic,jc) - grid_.c(ic-1,jc)) / dx;
            double fx = sigma * kappa(ic,jc) * dcdx;
            grid_.u(i,j) += dt * fx / rho;
        }
    }
    for (int j = 1; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            int ic = std::min(std::max(i,1), nx-1);
            int jc = std::min(std::max(j,1), ny-1);
            double dcdy = (grid_.c(ic,jc) - grid_.c(ic,jc-1)) / dy;
            double fy = sigma * kappa(ic,jc) * dcdy;
            grid_.v(i,j) += dt * fy / rho;
        }
    }
}

void Solver2D::step() {
    // CFL guard (simple): dt <= 0.9 * min(dx/|u|, dy/|v|)
    double umax = 1e-12, vmax = 1e-12;
    for (double v: grid_.u.data) umax = std::max(umax, std::abs(v));
    for (double v: grid_.v.data) vmax = std::max(vmax, std::abs(v));
    double cfl_dt = 0.9 * std::min(grid_.dx / umax, grid_.dy / vmax);
    if (params_.dt > cfl_dt) {
        // scale velocities down to satisfy CFL if user dt is too large (simple safeguard)
        double s1 = (umax>0)? (0.9*grid_.dx/(params_.dt*umax)) : 1.0;
        double s2 = (vmax>0)? (0.9*grid_.dy/(params_.dt*vmax)) : 1.0;
        double s = std::min(1.0, std::min(s1,s2));
        for (auto& x: grid_.u.data) x *= s;
        for (auto& x: grid_.v.data) x *= s;
    }

    advect_velocities();
    apply_body_forces();
    apply_surface_tension();
    apply_viscosity();
    project_incompressible();
    advect_vof_conservative();
    plic_reconstruct_and_correct();
}

double Solver2D::compute_mass() const {
    double sum = std::accumulate(grid_.c.data.begin(), grid_.c.data.end(), 0.0);
    return sum * grid_.dx * grid_.dy;
}

// Explicit template instantiations for common lambdas disabled; header-only templates suffice.

} // namespace vof2d
