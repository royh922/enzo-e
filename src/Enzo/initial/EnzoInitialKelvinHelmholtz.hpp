// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoInitialKelvinHelmholtz.hpp
/// @author   Generated for enzo-e based on Athena++ KH instability problem
/// @date     2025-09-12
/// @brief    [\ref Enzo] Initialization routine for Kelvin-Helmholtz instability
///           with cylindrical shear flow in hot diffuse medium

#ifndef ENZO_ENZO_INITIAL_KELVIN_HELMHOLTZ_HPP
#define ENZO_ENZO_INITIAL_KELVIN_HELMHOLTZ_HPP

class EnzoInitialKelvinHelmholtz : public Initial {
  /// @class    EnzoInitialKelvinHelmholtz
  /// @ingroup  Enzo
  /// @brief    [\ref Enzo] Initializer for Kelvin-Helmholtz instability
  /// 
  /// This class initializes a cylindrical cold gas region within hot diffuse
  /// medium with shear flow along the cylinder axis (x-direction), leading to
  /// Kelvin-Helmholtz instability. The setup includes:
  /// - Cylindrical cold gas region along x-axis with density contrast
  /// - Shear velocity profile with tanh transition
  /// - Optional velocity perturbations to seed the instability
  /// - Smoothing parameters for realistic transition regions

public: // interface

  /// Constructor
  EnzoInitialKelvinHelmholtz(int cycle, double time, ParameterGroup p) noexcept
    : Initial(cycle, time),
      gamma_adi_(p.value_float("gamma", 5.0/3.0)),
      n_0_(p.value_float("n_0", 0.01)),
      temperature_(p.value_float("temperature", 1.5)),
      density_contrast_(p.value_float("density_contrast", 10.0)),
      M_b_(p.value_float("M_b", 1.0)),
      radius_(p.value_float("radius", 0.5)),
      smoothing_thickness_(p.value_float("smoothing_thickness", 0.1)),
      smoothing_thickness_vel_(p.value_float("smoothing_thickness_vel", -1.0)),
      vel_pert_(p.value_float("vel_pert", 0.1)),
      lambda_pert_(p.value_float("lambda_pert", 1.0)),
      pert_loc_(p.value_float("pert_loc", 0.0)),
      uniform_bfield_{0.0, 0.0, 0.0},
      initialize_uniform_bfield_(false),
      noisy_ic_(p.value_logical("noisy_ic", false)),
      random_seed_(p.value_integer("random_seed", 12345))
  {
    // Use smoothing_thickness for velocity if not specified
    if (smoothing_thickness_vel_ < 0.0) {
      smoothing_thickness_vel_ = smoothing_thickness_;
    }

    // Handle uniform magnetic field initialization
    int uniform_bfield_length = p.list_length("uniform_bfield");
    if (uniform_bfield_length == 0) {
      initialize_uniform_bfield_ = false;
    } else if (uniform_bfield_length == 3) {
      initialize_uniform_bfield_ = true;
      for (int i = 0; i < 3; i++) {
        uniform_bfield_[i] = p.list_value_float(i, "uniform_bfield");
      }
    } else {
      ERROR("EnzoInitialKelvinHelmholtz",
            "Initial:kelvin_helmholtz:uniform_bfield must contain 0 or 3 entries.");
    }

    // Validate parameters
    ASSERT("EnzoInitialKelvinHelmholtz", "n_0 must be positive", n_0_ > 0);
    ASSERT("EnzoInitialKelvinHelmholtz", "temperature must be positive", temperature_ > 0);
    ASSERT("EnzoInitialKelvinHelmholtz", "density_contrast must be positive", 
           density_contrast_ > 0);
    ASSERT("EnzoInitialKelvinHelmholtz", "radius must be positive", radius_ > 0);
    ASSERT("EnzoInitialKelvinHelmholtz", "smoothing_thickness must be positive", 
           smoothing_thickness_ > 0);
    ASSERT("EnzoInitialKelvinHelmholtz", "smoothing_thickness_vel must be positive", 
           smoothing_thickness_vel_ > 0);
    ASSERT("EnzoInitialKelvinHelmholtz", "gamma must be greater than 1", 
           gamma_adi_ > 1.0);
  }

  /// CHARM++ PUP::able declaration
  PUPable_decl(EnzoInitialKelvinHelmholtz);

  /// CHARM++ migration constructor
  EnzoInitialKelvinHelmholtz(CkMigrateMessage *m)
    : Initial(m),
      gamma_adi_(5.0/3.0),
      n_0_(0.01),
      temperature_(1.5),
      density_contrast_(10.0),
      M_b_(1.0),
      radius_(0.5),
      smoothing_thickness_(0.1),
      smoothing_thickness_vel_(0.1),
      vel_pert_(0.1),
      lambda_pert_(1.0),
      pert_loc_(0.0),
      initialize_uniform_bfield_(false),
      noisy_ic_(false),
      random_seed_(12345)
  {
    for (int i = 0; i < 3; i++) { uniform_bfield_[i] = 0.0; }
  }

  /// CHARM++ Pack / Unpack function
  void pup(PUP::er &p)
  {
    // NOTE: update whenever attributes change
    TRACEPUP;

    Initial::pup(p);
    p | gamma_adi_;
    p | n_0_;
    p | temperature_;
    p | density_contrast_;
    p | M_b_;
    p | radius_;
    p | smoothing_thickness_;
    p | smoothing_thickness_vel_;
    p | vel_pert_;
    p | lambda_pert_;
    p | pert_loc_;
    p | initialize_uniform_bfield_;
    PUParray(p, uniform_bfield_, 3);
    p | noisy_ic_;
    p | random_seed_;
  }

public: // virtual methods

  /// Initialize the block
  virtual void enforce_block
  ( Block * block, const Hierarchy * hierarchy ) throw();

private: // attributes

  // NOTE: change pup() function whenever attributes change

  /// Adiabatic index (gamma)
  double gamma_adi_;

  /// Background number density (hot gas number density) [cm^-3]
  double n_0_;

  /// Stream temperature [10^4 K]
  double temperature_;

  /// Density contrast (cold gas density = n_0 * density_contrast)
  double density_contrast_;

  /// Stream Mach number relative to sound speed
  double M_b_;

  /// Radius of the cylindrical cold gas region
  double radius_;

  /// Smoothing thickness for density transition
  double smoothing_thickness_;

  /// Smoothing thickness for velocity perturbations
  double smoothing_thickness_vel_;

  /// Perturbation velocity amplitude
  double vel_pert_;

  /// Wavelength of perturbations (lambda_pert > 0 for sinusoidal, = 0 for localized)
  double lambda_pert_;

  /// Location of localized perturbation (used when lambda_pert = 0)
  double pert_loc_;

  /// Whether to initialize uniform magnetic field
  bool initialize_uniform_bfield_;

  /// Values of the uniform magnetic field
  double uniform_bfield_[3];

  /// Whether to add random noise to initial conditions
  bool noisy_ic_;

  /// Random seed for noise generation
  int random_seed_;
};

#endif // ENZO_ENZO_INITIAL_KELVIN_HELMHOLTZ_HPP