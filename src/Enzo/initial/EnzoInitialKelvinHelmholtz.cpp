// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialKelvinHelmholtz.cpp
/// @author   Generated for enzo-e based on Athena++ KH instability problem
/// @date     2025-09-12
/// @brief    [\ref Enzo] Implementation of EnzoInitialKelvinHelmholtz for
///           initializing Kelvin-Helmholtz instability with cylindrical shear flow

#include <cmath>
#include <random>

#include "Enzo/initial/initial.hpp"
#include "Enzo/enzo.hpp"
#include "Cello/cello.hpp"

//----------------------------------------------------------------------

void EnzoInitialKelvinHelmholtz::enforce_block
(Block * block, const Hierarchy * hierarchy) throw()
{
  Field field = block->data()->field();

  // Get field arrays
  EFlt3DArray density = field.view<enzo_float>("density");
  EFlt3DArray velocity_x = field.view<enzo_float>("velocity_x");
  EFlt3DArray velocity_y = field.view<enzo_float>("velocity_y");
  EFlt3DArray velocity_z = field.view<enzo_float>("velocity_z");
  EFlt3DArray total_energy = field.view<enzo_float>("total_energy");

  // Optional fields
  EFlt3DArray pressure;
  bool has_pressure = field.is_field("pressure");
  if (has_pressure) {
    pressure = field.view<enzo_float>("pressure");
  }

  // Handle magnetic fields
  EFlt3DArray bfield_x, bfield_y, bfield_z;
  bool has_bfield = (field.is_field("bfield_x") && 
                     field.is_field("bfield_y") && 
                     field.is_field("bfield_z"));
  if (has_bfield) {
    bfield_x = field.view<enzo_float>("bfield_x");
    bfield_y = field.view<enzo_float>("bfield_y");
    bfield_z = field.view<enzo_float>("bfield_z");
  }

  // Initialize random number generator for noise
  // Is this enough to ensure different noise on different blocks?
  std::minstd_rand generator(random_seed_);

  // Get grid dimensions and coordinates
  const int mx = density.shape(2);
  const int my = density.shape(1);
  const int mz = density.shape(0);

  // Get block bounds
  double xm, ym, zm;
  double xp, yp, zp;
  block->lower(&xm, &ym, &zm);
  block->upper(&xp, &yp, &zp);

  // Grid spacing
  const double hx = (xp - xm) / mx;
  const double hy = (yp - ym) / my;
  const double hz = (zp - zm) / mz;

  // Constants
  const double pi = cello::pi;
  const double rho_hot = rho_0_;
  const double rho_cold = rho_0_ * density_contrast_;

  // Initialize uniform magnetic field if requested
  if (has_bfield && initialize_uniform_bfield_) {
    for (int iz = 0; iz < mz; iz++) {
      for (int iy = 0; iy < my; iy++) {
        for (int ix = 0; ix < mx; ix++) {
          bfield_x(iz, iy, ix) = uniform_bfield_[0];
          bfield_y(iz, iy, ix) = uniform_bfield_[1];
          bfield_z(iz, iy, ix) = uniform_bfield_[2];
        }
      }
    }
  }

  // Calculate magnetic energy density if magnetic fields are present
  double magnetic_edens = 0.0;
  if (has_bfield) {
    magnetic_edens = 0.5 * (uniform_bfield_[0] * uniform_bfield_[0] + 
                           uniform_bfield_[1] * uniform_bfield_[1] + 
                           uniform_bfield_[2] * uniform_bfield_[2]);
  }

  // Main initialization loop
  for (int iz = 0; iz < mz; iz++) {
    // z coordinate at cell center
    double z = zm + (iz + 0.5) * hz;

    for (int iy = 0; iy < my; iy++) {
      // y coordinate at cell center  
      double y = ym + (iy + 0.5) * hy;

      for (int ix = 0; ix < mx; ix++) {
        // x coordinate at cell center
        double x = xm + (ix + 0.5) * hx;

        // Calculate radius from cylinder axis (x-axis)
        // Cylinder axis is along x, so radius is distance in y-z plane
        double r = std::sqrt(y * y + z * z);

        // Initialize density with tanh profile for smooth transition
        double density_val = rho_hot * (density_contrast_ / 2.0 + 0.5 + 
                            (density_contrast_ - 1.0) * 0.5 * 
                            (-std::tanh((r - radius_) / smoothing_thickness_)));

        density(iz, iy, ix) = density_val;

        // Initialize x-velocity (shear along cylinder axis)
        double vel_x = vel_shear_ * (-std::tanh((r - radius_) / smoothing_thickness_vel_));
        velocity_x(iz, iy, ix) = vel_x;

        // Initialize y and z velocities to zero initially
        double vel_y = 0.0;
        double vel_z = 0.0;

        // Add perturbations in the transition region
        if ((density_val > rho_hot) && (density_val < rho_cold)) {
          double mag = vel_pert_;
          
          // Apply Gaussian envelope centered on the cylinder boundary
          mag *= std::exp(-std::pow((r - radius_) / smoothing_thickness_vel_, 2));

          if (lambda_pert_ > 0.0) {
            // Sinusoidal perturbation along x-axis
            mag *= std::sin(2.0 * pi * x / lambda_pert_);
          } else if (lambda_pert_ == 0.0) {
            // Localized perturbation
            double pert_width = smoothing_thickness_vel_;
            mag *= std::exp(-std::pow((x - pert_loc_) / pert_width, 2));
          }

          // Calculate perturbation direction (cylindrical coordinates)
          if (r > 0.0) {
            double theta = std::atan2(y, z);
            vel_z += mag * std::cos(theta);  // z component
            vel_y += mag * std::sin(theta);  // y component
          }
        }

        // Add random noise if requested
        if (noisy_ic_) {
          std::uniform_real_distribution<double> dist(0.0, 1.0);
          vel_y *= dist(generator);
          vel_z *= dist(generator);
        }

        velocity_y(iz, iy, ix) = vel_y;
        velocity_z(iz, iy, ix) = vel_z;

        // Set pressure (constant throughout domain)
        if (has_pressure) {
          pressure(iz, iy, ix) = pgas_0_;
        }

        // Calculate total energy
        double kinetic_energy = 0.5 * (vel_x * vel_x + vel_y * vel_y + vel_z * vel_z);
        double thermal_energy = pgas_0_ / (gamma_adi_ - 1.0);
        
        double total_energy_val = thermal_energy + kinetic_energy;
        
        // Add magnetic energy if present
        if (has_bfield) {
          total_energy_val += magnetic_edens;
        }

        total_energy(iz, iy, ix) = total_energy_val;
      }
    }
  }

  block->initial_done();
}