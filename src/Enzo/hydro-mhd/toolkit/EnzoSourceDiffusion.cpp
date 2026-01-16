// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoSourceDiffusion.cpp
/// @author   Jui-Teng (Roy) Hsu (hjuiteng@gmail.com)
/// @brief    [\ref Enzo] Implementation of Enzo's EnzoSourceDiffusion class
///
/// Implements the full Navier-Stokes viscous stress tensor and thermal conduction
/// for compressible flows using a flux-based approach following Athena++.
///
/// The viscous stress tensor for an isotropic Newtonian fluid is:
///   tau_ij = mu * (dv_i/dx_j + dv_j/dx_i - (2/3) * delta_ij * div(v))
///
/// This is computed at cell faces and the divergence is used to update
/// conserved momentum and energy.
///
/// Performance notes:
/// - div(v) is precomputed at cell centers before the main loop
/// - Loop bounds are computed once outside loops
/// - EOS checks are done once and cached as booleans

#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"
#include "Enzo/hydro-mhd/hydro-mhd.hpp"

//----------------------------------------------------------------------

void EnzoSourceDiffusion::calculate_source(const double cur_dt, const EnzoEFltArrayMap& primitive_map,
                                           EnzoEFltArrayMap& dUcons_map,
                                           const std::array<enzo_float, 3>& cell_widths_xyz,
                                           const int stale_depth) const noexcept {
    // Only proceed if we have viscosity or conduction
    if ((viscosity_nu_ <= 0.0) && (thermal_kappa_ <= 0.0)) {
        return;
    }

    // Load primitive fields
    CelloView<const enzo_float, 3> rho = primitive_map.at("density");
    CelloView<const enzo_float, 3> vx = primitive_map.at("velocity_x");
    CelloView<const enzo_float, 3> vy = primitive_map.at("velocity_y");
    CelloView<const enzo_float, 3> vz = primitive_map.at("velocity_z");
    CelloView<const enzo_float, 3> p = primitive_map.at("pressure");

    // Accumulators for conserved quantity updates
    EFlt3DArray dU_vx = dUcons_map.at("velocity_x");
    EFlt3DArray dU_vy = dUcons_map.at("velocity_y");
    EFlt3DArray dU_vz = dUcons_map.at("velocity_z");
    EFlt3DArray dU_etot = dUcons_map.at("total_energy");

    // Check if we need to update internal energy (dual energy formalism)
    // Cache this check to avoid repeated virtual calls
    const bool use_dual_energy = enzo::fluid_props()->dual_energy_config().any_enabled();
    EFlt3DArray dU_eint;
    if (use_dual_energy) {
        dU_eint = dUcons_map.at("internal_energy");
    }

    // Array dimensions
    const int mz = rho.shape(0);
    const int my = rho.shape(1);
    const int mx = rho.shape(2);

    // Check which dimensions are active (cache these checks)
    const bool active_x = (mx > 1);
    const bool active_y = (my > 1);
    const bool active_z = (mz > 1);

    // Grid spacings - precompute inverse values
    const enzo_float dx = cell_widths_xyz[0];
    const enzo_float dy = cell_widths_xyz[1];
    const enzo_float dz = cell_widths_xyz[2];
    const enzo_float idx = active_x ? 1.0 / dx : 0.0;
    const enzo_float idy = active_y ? 1.0 / dy : 0.0;
    const enzo_float idz = active_z ? 1.0 / dz : 0.0;

    // Precompute commonly used derivative factors
    const enzo_float half_idx = 0.5 * idx;
    const enzo_float half_idy = 0.5 * idy;
    const enzo_float half_idz = 0.5 * idz;
    const enzo_float quarter_idx = 0.25 * idx;
    const enzo_float quarter_idy = 0.25 * idy;
    const enzo_float quarter_idz = 0.25 * idz;

    // EOS parameters - cache checks to avoid repeated virtual calls
    const bool ideal_eos = enzo::fluid_props()->eos_variant().holds_alternative<EnzoEOSIdeal>();
    enzo_float gm1 = 0.0;
    enzo_float inv_gm1 = 0.0;
    if (ideal_eos) {
        gm1 = enzo::fluid_props()->eos_variant().get<EnzoEOSIdeal>().gamma() - 1.0;
        inv_gm1 = 1.0 / gm1;
    }

    // Bulk viscosity coefficient: -2/3 for traceless stress tensor
    const enzo_float nuiso2 = -2.0 / 3.0;

    // Loop bounds for main computation:
    // We need to access neighbors at ix±1, iy±1, iz±1 for derivative computation,
    // so we need an extra layer of padding beyond stale_depth.
    // Following Enzo-E convention: loop from (stale_depth + 1) to (m - stale_depth - 1)
    // This ensures we don't access stale values in the ghost zones.
    const int iz_start = active_z ? stale_depth + 1 : 0;
    const int iz_end = active_z ? mz - stale_depth - 1 : 1;
    const int iy_start = active_y ? stale_depth + 1 : 0;
    const int iy_end = active_y ? my - stale_depth - 1 : 1;
    const int ix_start = active_x ? stale_depth + 1 : 0;
    const int ix_end = active_x ? mx - stale_depth - 1 : 1;
    //==========================================================================
    // PRECOMPUTE SCRATCH ARRAYS
    // Following Athena++, precompute div(v) at cell centers to avoid redundant
    // computation in the main loop. This is a significant optimization since
    // div(v) is needed at multiple face locations per cell.
    //==========================================================================

    // Allocate scratch array for div(v) - only if viscosity is enabled
    // The main loop runs from (stale_depth+1) to (m-stale_depth-1) and accesses
    // div_vel at ix±1, iy±1, iz±1 for face averaging. So div_vel must be valid
    // from stale_depth to (m-stale_depth-1).
    // 
    // div_vel computation uses centered differences requiring ±1 neighbors,
    // BUT the primitive velocities (vx, vy, vz) are still valid in ghost zones
    // since source terms don't increase stale_depth (only flux divergence does).
    // So we can safely compute div_vel from stale_depth to (m-stale_depth).
    EFlt3DArray div_vel_arr;
    if (viscosity_nu_ > 0.0) {
        // Compute div_vel at range [stale_depth, m-stale_depth) to cover
        // all cells accessed by main loop's face-averaging operations
        const int div_iz_start = active_z ? stale_depth : 0;
        const int div_iz_end = active_z ? mz - stale_depth : 1;
        const int div_iy_start = active_y ? stale_depth : 0;
        const int div_iy_end = active_y ? my - stale_depth : 1;
        const int div_ix_start = active_x ? stale_depth : 0;
        const int div_ix_end = active_x ? mx - stale_depth : 1;
        div_vel_arr = EFlt3DArray(mz, my, mx);

        // Compute div(v) at all needed cell centers
        for (int iz = div_iz_start; iz < div_iz_end; iz++) {
            for (int iy = div_iy_start; iy < div_iy_end; iy++) {
                for (int ix = div_ix_start; ix < div_ix_end; ix++) {
                    enzo_float div = 0.0;
                    if (active_x) {
                        div += (vx(iz, iy, ix + 1) - vx(iz, iy, ix - 1)) * half_idx;
                    }
                    if (active_y) {
                        div += (vy(iz, iy + 1, ix) - vy(iz, iy - 1, ix)) * half_idy;
                    }
                    if (active_z) {
                        div += (vz(iz + 1, iy, ix) - vz(iz - 1, iy, ix)) * half_idz;
                    }
                    div_vel_arr(iz, iy, ix) = div;
                }
            }
        }
    }

    //==========================================================================
    // 1. VISCOSITY: Full Navier-Stokes stress tensor
    //    tau_ij = mu * (dv_i/dx_j + dv_j/dx_i - (2/3)*delta_ij*div(v))
    //
    //    Following Athena++, we compute viscous fluxes at cell faces and
    //    use their divergence to update conserved momentum and energy.
    //==========================================================================

    if (viscosity_nu_ > 0.0) {
        for (int iz = iz_start; iz < iz_end; iz++) {
            for (int iy = iy_start; iy < iy_end; iy++) {
                for (int ix = ix_start; ix < ix_end; ix++) {
                    // Get precomputed divergence at this cell
                    const enzo_float div_c = div_vel_arr(iz, iy, ix);

                    // Dynamic viscosity at cell center
                    const enzo_float mu_c = rho(iz, iy, ix) * viscosity_nu_;

                    // Accumulate momentum and energy source terms
                    enzo_float src_momx = 0.0;
                    enzo_float src_momy = 0.0;
                    enzo_float src_momz = 0.0;
                    enzo_float src_etot = 0.0;

                    //------------------------------------------------------------------
                    // X-direction fluxes (at i+1/2 and i-1/2 faces)
                    //------------------------------------------------------------------
                    if (active_x) {
                        // Right face (i+1/2)
                        const enzo_float mu_r = 0.5 * (mu_c + rho(iz, iy, ix + 1) * viscosity_nu_);
                        const enzo_float div_r = 0.5 * (div_c + div_vel_arr(iz, iy, ix + 1));

                        // Velocity derivatives at right face
                        const enzo_float dvx_dx_r = (vx(iz, iy, ix + 1) - vx(iz, iy, ix)) * idx;
                        const enzo_float dvy_dx_r = (vy(iz, iy, ix + 1) - vy(iz, iy, ix)) * idx;
                        const enzo_float dvz_dx_r = (vz(iz, iy, ix + 1) - vz(iz, iy, ix)) * idx;

                        enzo_float dvx_dy_r = 0.0;
                        if (active_y) {
                            dvx_dy_r = ((vx(iz, iy + 1, ix + 1) - vx(iz, iy - 1, ix + 1)) +
                                        (vx(iz, iy + 1, ix) - vx(iz, iy - 1, ix))) *
                                       quarter_idy;
                        }
                        enzo_float dvx_dz_r = 0.0;
                        if (active_z) {
                            dvx_dz_r = ((vx(iz + 1, iy, ix + 1) - vx(iz - 1, iy, ix + 1)) +
                                        (vx(iz + 1, iy, ix) - vx(iz - 1, iy, ix))) *
                                       quarter_idz;
                        }

                        // Viscous stress components at right face
                        const enzo_float tau_xx_r = mu_r * (2.0 * dvx_dx_r + nuiso2 * div_r);
                        const enzo_float tau_xy_r = mu_r * (dvx_dy_r + dvy_dx_r);
                        const enzo_float tau_xz_r = mu_r * (dvx_dz_r + dvz_dx_r);

                        // Left face (i-1/2)
                        const enzo_float mu_l = 0.5 * (mu_c + rho(iz, iy, ix - 1) * viscosity_nu_);
                        const enzo_float div_l = 0.5 * (div_c + div_vel_arr(iz, iy, ix - 1));

                        const enzo_float dvx_dx_l = (vx(iz, iy, ix) - vx(iz, iy, ix - 1)) * idx;
                        const enzo_float dvy_dx_l = (vy(iz, iy, ix) - vy(iz, iy, ix - 1)) * idx;
                        const enzo_float dvz_dx_l = (vz(iz, iy, ix) - vz(iz, iy, ix - 1)) * idx;

                        enzo_float dvx_dy_l = 0.0;
                        if (active_y) {
                            dvx_dy_l = ((vx(iz, iy + 1, ix) - vx(iz, iy - 1, ix)) +
                                        (vx(iz, iy + 1, ix - 1) - vx(iz, iy - 1, ix - 1))) *
                                       quarter_idy;
                        }
                        enzo_float dvx_dz_l = 0.0;
                        if (active_z) {
                            dvx_dz_l = ((vx(iz + 1, iy, ix) - vx(iz - 1, iy, ix)) +
                                        (vx(iz + 1, iy, ix - 1) - vx(iz - 1, iy, ix - 1))) *
                                       quarter_idz;
                        }

                        const enzo_float tau_xx_l = mu_l * (2.0 * dvx_dx_l + nuiso2 * div_l);
                        const enzo_float tau_xy_l = mu_l * (dvx_dy_l + dvy_dx_l);
                        const enzo_float tau_xz_l = mu_l * (dvx_dz_l + dvz_dx_l);

                        // Flux divergence: d(tau)/dx (note: flux = -tau, so div gives +tau)
                        src_momx += (tau_xx_r - tau_xx_l) * idx;
                        src_momy += (tau_xy_r - tau_xy_l) * idx;
                        src_momz += (tau_xz_r - tau_xz_l) * idx;

                        // Energy flux: v . tau at face (work done by viscous stress)
                        if (ideal_eos) {
                            const enzo_float vx_r = 0.5 * (vx(iz, iy, ix) + vx(iz, iy, ix + 1));
                            const enzo_float vy_r = 0.5 * (vy(iz, iy, ix) + vy(iz, iy, ix + 1));
                            const enzo_float vz_r = 0.5 * (vz(iz, iy, ix) + vz(iz, iy, ix + 1));
                            const enzo_float en_flux_r = vx_r * tau_xx_r + vy_r * tau_xy_r + vz_r * tau_xz_r;

                            const enzo_float vx_l = 0.5 * (vx(iz, iy, ix) + vx(iz, iy, ix - 1));
                            const enzo_float vy_l = 0.5 * (vy(iz, iy, ix) + vy(iz, iy, ix - 1));
                            const enzo_float vz_l = 0.5 * (vz(iz, iy, ix) + vz(iz, iy, ix - 1));
                            const enzo_float en_flux_l = vx_l * tau_xx_l + vy_l * tau_xy_l + vz_l * tau_xz_l;

                            src_etot += (en_flux_r - en_flux_l) * idx;
                        }
                    }

                    //------------------------------------------------------------------
                    // Y-direction fluxes (at j+1/2 and j-1/2 faces)
                    //------------------------------------------------------------------
                    if (active_y) {
                        // Top face (j+1/2)
                        const enzo_float mu_t = 0.5 * (mu_c + rho(iz, iy + 1, ix) * viscosity_nu_);
                        const enzo_float div_t = 0.5 * (div_c + div_vel_arr(iz, iy + 1, ix));

                        const enzo_float dvy_dy_t = (vy(iz, iy + 1, ix) - vy(iz, iy, ix)) * idy;
                        const enzo_float dvx_dy_t = (vx(iz, iy + 1, ix) - vx(iz, iy, ix)) * idy;
                        const enzo_float dvz_dy_t = (vz(iz, iy + 1, ix) - vz(iz, iy, ix)) * idy;

                        enzo_float dvy_dx_t = 0.0;
                        if (active_x) {
                            dvy_dx_t = ((vy(iz, iy + 1, ix + 1) - vy(iz, iy + 1, ix - 1)) +
                                        (vy(iz, iy, ix + 1) - vy(iz, iy, ix - 1))) *
                                       quarter_idx;
                        }
                        enzo_float dvy_dz_t = 0.0;
                        if (active_z) {
                            dvy_dz_t = ((vy(iz + 1, iy + 1, ix) - vy(iz - 1, iy + 1, ix)) +
                                        (vy(iz + 1, iy, ix) - vy(iz - 1, iy, ix))) *
                                       quarter_idz;
                        }

                        const enzo_float tau_yy_t = mu_t * (2.0 * dvy_dy_t + nuiso2 * div_t);
                        const enzo_float tau_yx_t = mu_t * (dvy_dx_t + dvx_dy_t);
                        const enzo_float tau_yz_t = mu_t * (dvy_dz_t + dvz_dy_t);

                        // Bottom face (j-1/2)
                        const enzo_float mu_b = 0.5 * (mu_c + rho(iz, iy - 1, ix) * viscosity_nu_);
                        const enzo_float div_b = 0.5 * (div_c + div_vel_arr(iz, iy - 1, ix));

                        const enzo_float dvy_dy_b = (vy(iz, iy, ix) - vy(iz, iy - 1, ix)) * idy;
                        const enzo_float dvx_dy_b = (vx(iz, iy, ix) - vx(iz, iy - 1, ix)) * idy;
                        const enzo_float dvz_dy_b = (vz(iz, iy, ix) - vz(iz, iy - 1, ix)) * idy;

                        enzo_float dvy_dx_b = 0.0;
                        if (active_x) {
                            dvy_dx_b = ((vy(iz, iy, ix + 1) - vy(iz, iy, ix - 1)) +
                                        (vy(iz, iy - 1, ix + 1) - vy(iz, iy - 1, ix - 1))) *
                                       quarter_idx;
                        }
                        enzo_float dvy_dz_b = 0.0;
                        if (active_z) {
                            dvy_dz_b = ((vy(iz + 1, iy, ix) - vy(iz - 1, iy, ix)) +
                                        (vy(iz + 1, iy - 1, ix) - vy(iz - 1, iy - 1, ix))) *
                                       quarter_idz;
                        }

                        const enzo_float tau_yy_b = mu_b * (2.0 * dvy_dy_b + nuiso2 * div_b);
                        const enzo_float tau_yx_b = mu_b * (dvy_dx_b + dvx_dy_b);
                        const enzo_float tau_yz_b = mu_b * (dvy_dz_b + dvz_dy_b);

                        src_momx += (tau_yx_t - tau_yx_b) * idy;
                        src_momy += (tau_yy_t - tau_yy_b) * idy;
                        src_momz += (tau_yz_t - tau_yz_b) * idy;

                        if (ideal_eos) {
                            const enzo_float vx_t = 0.5 * (vx(iz, iy, ix) + vx(iz, iy + 1, ix));
                            const enzo_float vy_t = 0.5 * (vy(iz, iy, ix) + vy(iz, iy + 1, ix));
                            const enzo_float vz_t = 0.5 * (vz(iz, iy, ix) + vz(iz, iy + 1, ix));
                            const enzo_float en_flux_t = vx_t * tau_yx_t + vy_t * tau_yy_t + vz_t * tau_yz_t;

                            const enzo_float vx_b = 0.5 * (vx(iz, iy, ix) + vx(iz, iy - 1, ix));
                            const enzo_float vy_b = 0.5 * (vy(iz, iy, ix) + vy(iz, iy - 1, ix));
                            const enzo_float vz_b = 0.5 * (vz(iz, iy, ix) + vz(iz, iy - 1, ix));
                            const enzo_float en_flux_b = vx_b * tau_yx_b + vy_b * tau_yy_b + vz_b * tau_yz_b;

                            src_etot += (en_flux_t - en_flux_b) * idy;
                        }
                    }

                    //------------------------------------------------------------------
                    // Z-direction fluxes (at k+1/2 and k-1/2 faces)
                    //------------------------------------------------------------------
                    if (active_z) {
                        // Front face (k+1/2)
                        const enzo_float mu_f = 0.5 * (mu_c + rho(iz + 1, iy, ix) * viscosity_nu_);
                        const enzo_float div_f = 0.5 * (div_c + div_vel_arr(iz + 1, iy, ix));

                        const enzo_float dvz_dz_f = (vz(iz + 1, iy, ix) - vz(iz, iy, ix)) * idz;
                        const enzo_float dvx_dz_f = (vx(iz + 1, iy, ix) - vx(iz, iy, ix)) * idz;
                        const enzo_float dvy_dz_f = (vy(iz + 1, iy, ix) - vy(iz, iy, ix)) * idz;

                        enzo_float dvz_dx_f = 0.0;
                        if (active_x) {
                            dvz_dx_f = ((vz(iz + 1, iy, ix + 1) - vz(iz + 1, iy, ix - 1)) +
                                        (vz(iz, iy, ix + 1) - vz(iz, iy, ix - 1))) *
                                       quarter_idx;
                        }
                        enzo_float dvz_dy_f = 0.0;
                        if (active_y) {
                            dvz_dy_f = ((vz(iz + 1, iy + 1, ix) - vz(iz + 1, iy - 1, ix)) +
                                        (vz(iz, iy + 1, ix) - vz(iz, iy - 1, ix))) *
                                       quarter_idy;
                        }

                        const enzo_float tau_zz_f = mu_f * (2.0 * dvz_dz_f + nuiso2 * div_f);
                        const enzo_float tau_zx_f = mu_f * (dvz_dx_f + dvx_dz_f);
                        const enzo_float tau_zy_f = mu_f * (dvz_dy_f + dvy_dz_f);

                        // Back face (k-1/2)
                        const enzo_float mu_bk = 0.5 * (mu_c + rho(iz - 1, iy, ix) * viscosity_nu_);
                        const enzo_float div_bk = 0.5 * (div_c + div_vel_arr(iz - 1, iy, ix));

                        const enzo_float dvz_dz_bk = (vz(iz, iy, ix) - vz(iz - 1, iy, ix)) * idz;
                        const enzo_float dvx_dz_bk = (vx(iz, iy, ix) - vx(iz - 1, iy, ix)) * idz;
                        const enzo_float dvy_dz_bk = (vy(iz, iy, ix) - vy(iz - 1, iy, ix)) * idz;

                        enzo_float dvz_dx_bk = 0.0;
                        if (active_x) {
                            dvz_dx_bk = ((vz(iz, iy, ix + 1) - vz(iz, iy, ix - 1)) +
                                         (vz(iz - 1, iy, ix + 1) - vz(iz - 1, iy, ix - 1))) *
                                        quarter_idx;
                        }
                        enzo_float dvz_dy_bk = 0.0;
                        if (active_y) {
                            dvz_dy_bk = ((vz(iz, iy + 1, ix) - vz(iz, iy - 1, ix)) +
                                         (vz(iz - 1, iy + 1, ix) - vz(iz - 1, iy - 1, ix))) *
                                        quarter_idy;
                        }

                        const enzo_float tau_zz_bk = mu_bk * (2.0 * dvz_dz_bk + nuiso2 * div_bk);
                        const enzo_float tau_zx_bk = mu_bk * (dvz_dx_bk + dvx_dz_bk);
                        const enzo_float tau_zy_bk = mu_bk * (dvz_dy_bk + dvy_dz_bk);

                        src_momx += (tau_zx_f - tau_zx_bk) * idz;
                        src_momy += (tau_zy_f - tau_zy_bk) * idz;
                        src_momz += (tau_zz_f - tau_zz_bk) * idz;

                        if (ideal_eos) {
                            const enzo_float vx_f = 0.5 * (vx(iz, iy, ix) + vx(iz + 1, iy, ix));
                            const enzo_float vy_f = 0.5 * (vy(iz, iy, ix) + vy(iz + 1, iy, ix));
                            const enzo_float vz_f = 0.5 * (vz(iz, iy, ix) + vz(iz + 1, iy, ix));
                            const enzo_float en_flux_f = vx_f * tau_zx_f + vy_f * tau_zy_f + vz_f * tau_zz_f;

                            const enzo_float vx_bk = 0.5 * (vx(iz, iy, ix) + vx(iz - 1, iy, ix));
                            const enzo_float vy_bk = 0.5 * (vy(iz, iy, ix) + vy(iz - 1, iy, ix));
                            const enzo_float vz_bk = 0.5 * (vz(iz, iy, ix) + vz(iz - 1, iy, ix));
                            const enzo_float en_flux_bk = vx_bk * tau_zx_bk + vy_bk * tau_zy_bk + vz_bk * tau_zz_bk;

                            src_etot += (en_flux_f - en_flux_bk) * idz;
                        }
                    }

                    //------------------------------------------------------------------
                    // Viscous dissipation for internal energy (dual energy formalism)
                    // Phi = tau_ij * S_ij (always positive, heats the gas)
                    //------------------------------------------------------------------
                    enzo_float src_eint = 0.0;
                    if (ideal_eos && use_dual_energy) {
                        // Compute strain rate tensor components at cell center
                        enzo_float dvx_dx = 0.0, dvy_dy = 0.0, dvz_dz = 0.0;
                        enzo_float dvx_dy = 0.0, dvx_dz = 0.0;
                        enzo_float dvy_dx = 0.0, dvy_dz = 0.0;
                        enzo_float dvz_dx = 0.0, dvz_dy = 0.0;

                        if (active_x) {
                            dvx_dx = (vx(iz, iy, ix + 1) - vx(iz, iy, ix - 1)) * half_idx;
                            dvy_dx = (vy(iz, iy, ix + 1) - vy(iz, iy, ix - 1)) * half_idx;
                            dvz_dx = (vz(iz, iy, ix + 1) - vz(iz, iy, ix - 1)) * half_idx;
                        }
                        if (active_y) {
                            dvx_dy = (vx(iz, iy + 1, ix) - vx(iz, iy - 1, ix)) * half_idy;
                            dvy_dy = (vy(iz, iy + 1, ix) - vy(iz, iy - 1, ix)) * half_idy;
                            dvz_dy = (vz(iz, iy + 1, ix) - vz(iz, iy - 1, ix)) * half_idy;
                        }
                        if (active_z) {
                            dvx_dz = (vx(iz + 1, iy, ix) - vx(iz - 1, iy, ix)) * half_idz;
                            dvy_dz = (vy(iz + 1, iy, ix) - vy(iz - 1, iy, ix)) * half_idz;
                            dvz_dz = (vz(iz + 1, iy, ix) - vz(iz - 1, iy, ix)) * half_idz;
                        }

                        // Strain rate tensor: S_ij = 0.5 * (dv_i/dx_j + dv_j/dx_i)
                        const enzo_float S_xx = dvx_dx;
                        const enzo_float S_yy = dvy_dy;
                        const enzo_float S_zz = dvz_dz;
                        const enzo_float S_xy = 0.5 * (dvx_dy + dvy_dx);
                        const enzo_float S_xz = 0.5 * (dvx_dz + dvz_dx);
                        const enzo_float S_yz = 0.5 * (dvy_dz + dvz_dy);

                        // Viscous dissipation: Phi = tau_ij * S_ij
                        // = 2*mu*(S_ij*S_ij) - (2/3)*mu*div(v)^2
                        const enzo_float S_sq =
                            S_xx * S_xx + S_yy * S_yy + S_zz * S_zz + 2.0 * (S_xy * S_xy + S_xz * S_xz + S_yz * S_yz);
                        src_eint = 2.0 * mu_c * S_sq + nuiso2 * mu_c * div_c * div_c;
                    }

                    // Apply updates
                    dU_vx(iz, iy, ix) += cur_dt * src_momx;
                    dU_vy(iz, iy, ix) += cur_dt * src_momy;
                    dU_vz(iz, iy, ix) += cur_dt * src_momz;
                    dU_etot(iz, iy, ix) += cur_dt * src_etot;
                    if (use_dual_energy) {
                        dU_eint(iz, iy, ix) += cur_dt * src_eint;
                    }

                }  // ix
            }  // iy
        }  // iz
    }  // viscosity

    //==========================================================================
    // 2. THERMAL CONDUCTION (Fourier's Law)
    //    Heat flux: q = -kappa * grad(T) = -rho * chi * grad(e)
    //    d/dt(rho*e) = -div(q) = div(rho * chi * grad(e))
    //==========================================================================

    if ((thermal_kappa_ > 0.0) && ideal_eos) {
        for (int iz = iz_start; iz < iz_end; iz++) {
            for (int iy = iy_start; iy < iy_end; iy++) {
                for (int ix = ix_start; ix < ix_end; ix++) {
                    // Specific internal energy at cell center: e = p / (rho * (gamma-1))
                    const enzo_float rho_c = rho(iz, iy, ix);
                    const enzo_float e_c = p(iz, iy, ix) * inv_gm1 / rho_c;
                    const enzo_float alpha_c = rho_c * thermal_kappa_;

                    enzo_float src_energy = 0.0;

                    // X-direction heat flux
                    if (active_x) {
                        const enzo_float rho_r = rho(iz, iy, ix + 1);
                        const enzo_float rho_l = rho(iz, iy, ix - 1);
                        const enzo_float alpha_r = 0.5 * (alpha_c + rho_r * thermal_kappa_);
                        const enzo_float alpha_l = 0.5 * (alpha_c + rho_l * thermal_kappa_);
                        const enzo_float e_r = p(iz, iy, ix + 1) * inv_gm1 / rho_r;
                        const enzo_float e_l = p(iz, iy, ix - 1) * inv_gm1 / rho_l;
                        const enzo_float de_dx_r = (e_r - e_c) * idx;
                        const enzo_float de_dx_l = (e_c - e_l) * idx;
                        src_energy += (alpha_r * de_dx_r - alpha_l * de_dx_l) * idx;
                    }

                    // Y-direction heat flux
                    if (active_y) {
                        const enzo_float rho_t = rho(iz, iy + 1, ix);
                        const enzo_float rho_b = rho(iz, iy - 1, ix);
                        const enzo_float alpha_t = 0.5 * (alpha_c + rho_t * thermal_kappa_);
                        const enzo_float alpha_b = 0.5 * (alpha_c + rho_b * thermal_kappa_);
                        const enzo_float e_t = p(iz, iy + 1, ix) * inv_gm1 / rho_t;
                        const enzo_float e_b = p(iz, iy - 1, ix) * inv_gm1 / rho_b;
                        const enzo_float de_dy_t = (e_t - e_c) * idy;
                        const enzo_float de_dy_b = (e_c - e_b) * idy;
                        src_energy += (alpha_t * de_dy_t - alpha_b * de_dy_b) * idy;
                    }

                    // Z-direction heat flux
                    if (active_z) {
                        const enzo_float rho_f = rho(iz + 1, iy, ix);
                        const enzo_float rho_bk = rho(iz - 1, iy, ix);
                        const enzo_float alpha_f = 0.5 * (alpha_c + rho_f * thermal_kappa_);
                        const enzo_float alpha_bk = 0.5 * (alpha_c + rho_bk * thermal_kappa_);
                        const enzo_float e_f = p(iz + 1, iy, ix) * inv_gm1 / rho_f;
                        const enzo_float e_bk = p(iz - 1, iy, ix) * inv_gm1 / rho_bk;
                        const enzo_float de_dz_f = (e_f - e_c) * idz;
                        const enzo_float de_dz_bk = (e_c - e_bk) * idz;
                        src_energy += (alpha_f * de_dz_f - alpha_bk * de_dz_bk) * idz;
                    }

                    dU_etot(iz, iy, ix) += cur_dt * src_energy;
                    if (use_dual_energy) {
                        dU_eint(iz, iy, ix) += cur_dt * src_energy;
                    }

                }  // ix
            }  // iy
        }  // iz
    }  // thermal conduction
}
