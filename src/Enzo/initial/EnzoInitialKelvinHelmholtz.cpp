// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialKelvinHelmholtz.cpp
/// @author   Generated for enzo-e based on Athena++ KH instability problem
/// @date     2025-09-12
/// @brief    [\ref Enzo] Implementation of EnzoInitialKelvinHelmholtz for
///           initializing Kelvin-Helmholtz instability with cylindrical shear flow
///
/// UNIT SYSTEM:
/// - Mass unit: Hydrogen atom mass (m_H = 1.673e-24 g)
/// - Length unit: kpc (3.086e21 cm)
/// - Time unit: derived from length and velocity units
/// - Temperature unit: 10^4 K
/// - Number density: cm^-3
/// - Mass density: number density * molecular weight (in code units)

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

  EFlt3DArray internal_energy;
  bool has_internal_energy = field.is_field("internal_energy");
  if (has_internal_energy) {
    internal_energy = field.view<enzo_float>("internal_energy");
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

  // Constants and unit conversions
  const double pi = cello::pi;
  const double m_H = enzo_constants::mass_hydrogen;  // Hydrogen mass in grams 
  const double k_B = enzo_constants::kboltz;     // Boltzmann constant [erg/K]

  // Read code unit conversions
  const EnzoUnits * enzo_units = enzo::units();
  
  // Get molecular weight from Enzo config
  const EnzoConfig * enzo_config = enzo::config();

  const double mol_weight = enzo_config->physics_fluid_props_mol_weight;
  
  
  // Calculate densities in code units (mass density = number density * molecular weight)
  // Stream = cold, dense gas in the center (inside cylinder)
  // Surrounding = hot, diffuse gas outside cylinder
  // Convert to code units
  const double rho_cold = n_0_ * mol_weight * m_H * std::pow(enzo_units->length(), 3) / enzo_units->mass();
  const double rho_hot = rho_cold / density_contrast_;  // Hot (surrounding) gas mass density
  
  // For uniform pressure: P = constant
  // From ideal gas law: P = ρ * c_s^2 / γ = ρ * k_B * T / (μ * m_H * γ)
  // With uniform pressure, higher density regions have lower temperature
  // This maintains pressure equilibrium across the interface
  const double T_cold = temperature_;  // Cold stream temperature in K
  // Uniform pressure in code units
  const double p_uniform = rho_cold * k_B * T_cold / (mol_weight * m_H) / std::pow(enzo_units->length() / enzo_units->time(), 2);  
  
  // Calculate temperatures for uniform pressure
  const double T_hot = T_cold * density_contrast_;  // Hot gas temperature
  // T_cold is already defined above
  
  // Calculate sound speed and shear velocity using cold stream properties
  const double c_b = std::sqrt(gamma_adi_ * p_uniform / rho_hot);  // Background sound speed
  const double vel_shear = M_b_ * c_b;  // Shear velocity based on Mach number

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

  // Calculate specific magnetic energy density if magnetic fields are present
  double magnetic_edens = 0.0;
  if (has_bfield) {
    magnetic_edens = 0.5 * (uniform_bfield_[0] * uniform_bfield_[0] + 
                           uniform_bfield_[1] * uniform_bfield_[1] + 
                           uniform_bfield_[2] * uniform_bfield_[2]);
  }

  // Initialize random number generator if needed
  std::mt19937 generator;
  if (noisy_ic_) {
    generator.seed(random_seed_);
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

        // Initialize density based on radius
        // Inside cylinder (r < radius): cold, dense stream
        // Outside cylinder (r >= radius): hot, diffuse surrounding gas
        double density_val = r < radius_ ? rho_cold : rho_hot;
        density(iz, iy, ix) = density_val;

        // Initialize x-velocity (shear along cylinder axis)
        // Cold stream moves with shear velocity, hot gas is at rest
        double vel_x = (r < radius_) ? vel_shear : 0.0;
        velocity_x(iz, iy, ix) = vel_x;

        // Initialize y and z velocities to zero initially
        double vel_y = 0.0;
        double vel_z = 0.0;

        // Should be multimode but for now just single mode
        // Add perturbations at the interface
        if (std::abs(r - radius_) < 0.5) {  // Near the interface
          double mag = vel_pert_ * c_b;
          
          if (lambda_pert_ > 0.0) {
            // Sinusoidal perturbation along x-axis
            mag *= std::sin(2.0 * pi * x / lambda_pert_);
          } else if (lambda_pert_ == 0.0) {
            // Localized perturbation
            double pert_width = 1.0;  // Can be made configurable
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
          std::uniform_real_distribution<double> dist(-0.1, 0.1);
          vel_y += dist(generator) * vel_pert_;
          vel_z += dist(generator) * vel_pert_;
        }

        velocity_y(iz, iy, ix) = vel_y;
        velocity_z(iz, iy, ix) = vel_z;

        // Set pressure (uniform throughout domain)
        if (has_pressure) {
          pressure(iz, iy, ix) = p_uniform;
        }

        // Note that enzo-e uses specific energy (per mass) in dealing with energy
        // Set internal energy (uniform pressure, but varies with density)
        if (has_internal_energy) {
          internal_energy(iz, iy, ix) = p_uniform / (gamma_adi_ - 1.0) / density_val;
        }

        // Calculate total energy
        double kinetic_energy = 0.5 * (vel_x * vel_x + vel_y * vel_y + vel_z * vel_z);
        double thermal_energy = p_uniform / (gamma_adi_ - 1.0) / density_val;
        
        double total_energy_val = thermal_energy + kinetic_energy;
        
        // Add specific magnetic energy if present
        if (has_bfield) {
          total_energy_val += magnetic_edens / density_val;
        }

        total_energy(iz, iy, ix) = total_energy_val;
      }
    }
  }

  block->initial_done();
}