// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSourceDiffusion.hpp
/// @author   Jui-Teng (Roy) Hsu (hjuiteng@gmail.com)
/// @date     Thu Jan 15 2026
/// @brief    [\ref Enzo] Declaration of Enzo's SourceDiffusion class. This
/// computes the momentum and energy source terms from physical diffusion
/// (viscosity and thermal conduction) for compressible flows.


#ifndef ENZO_ENZO_SOURCE_DIFFUSION_HPP
#define ENZO_ENZO_SOURCE_DIFFUSION_HPP

class EnzoSourceDiffusion
{
  /// @class    EnzoSourceDiffusion
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Computes diffusion source terms for hydro solvers
  ///
  /// Implements the full Navier-Stokes viscous stress tensor for compressible
  /// flows and Fourier thermal conduction, following the approach in Athena++.
  ///
  /// The viscous stress tensor for an isotropic Newtonian fluid is:
  ///   tau_ij = mu * (dv_i/dx_j + dv_j/dx_i - (2/3)*delta_ij*div(v))
  /// where mu = rho * nu is the dynamic viscosity.
  ///
  /// Viscous fluxes are computed at cell faces and their divergence updates
  /// the conserved momentum. The energy equation includes the work done by
  /// viscous stresses (v . div(tau)) and viscous dissipation (tau : grad(v)).

public:

  /// Constructor
  /// @param[in] viscosity_nu Kinematic viscosity coefficient [L^2/T]
  /// @param[in] thermal_kappa Thermal diffusivity coefficient [L^2/T]
  EnzoSourceDiffusion(enzo_float viscosity_nu, enzo_float thermal_kappa) noexcept
    : viscosity_nu_(viscosity_nu), thermal_kappa_(thermal_kappa)
  {}

  /// Computes the momentum and energy source terms from the full Navier-Stokes
  /// viscous stress tensor and thermal conduction (Fourier's law).
  ///
  /// The viscous stress tensor includes:
  ///   - Shear stress: mu * (dv_i/dx_j + dv_j/dx_i)
  ///   - Bulk viscosity: -(2/3) * mu * div(v) * delta_ij
  ///
  /// Source terms are added to arrays with keys "velocity_x", "velocity_y",
  /// "velocity_z", and "total_energy" (and "internal_energy" if dual energy
  /// formalism is enabled).
  ///
  /// @param[in]  cur_dt The timestep over which to apply the source term
  /// @param[in]  primitive_map Map holding cell-centered density, velocity,
  ///     and pressure.
  /// @param[out] dUcons_map Map of arrays where source term contributions to
  ///     conserved quantities are accumulated.
  /// @param[in]  cell_widths_xyz Array with cell widths [dx, dy, dz]
  /// @param[in]  stale_depth indicates the current stale_depth for cell-centered
  ///     quantities. Only cells with distance >= stale_depth from boundaries are
  ///     updated (to avoid accessing invalid ghost zones).
  void calculate_source(const double cur_dt,
                        const EnzoEFltArrayMap &primitive_map,
                        EnzoEFltArrayMap &dUcons_map,
                        const std::array<enzo_float,3> &cell_widths_xyz,
			const int stale_depth) const noexcept;

private:
  enzo_float viscosity_nu_;   ///< Kinematic viscosity [L^2/T]
  enzo_float thermal_kappa_;  ///< Thermal diffusivity [L^2/T]
};

#endif /* ENZO_ENZO_SOURCE_DIFFUSION_HPP */
