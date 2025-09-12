# Kelvin-Helmholtz Instability Initial Condition for Enzo-E

This document describes the implementation of the Kelvin-Helmholtz instability initial condition for Enzo-E, based on the Athena++ reference implementation.

## Overview

The Kelvin-Helmholtz instability is a fundamental fluid instability that occurs at the interface between two fluids with different velocities. This implementation creates a cylindrical cold gas region within a hot diffuse medium with shear flow along the cylinder axis, leading to the development of KH instabilities.

## Physical Setup

### Geometry
- **Cylindrical cold gas region**: Dense gas cylinder along the x-axis (cylinder axis)
- **Hot diffuse medium**: Lower density gas surrounding the cylinder
- **Coordinate system**: 
  - x-axis: Along the cylinder (shear direction)
  - y-z plane: Perpendicular to cylinder (defines radius r = √(y² + z²))

### Density Profile
The density is set using a hyperbolic tangent (tanh) transition:
```
ρ(r) = ρ₀ * [(ρ_contrast/2) + 0.5 + (ρ_contrast - 1.0) * 0.5 * (-tanh((r - R)/δ))]
```
Where:
- ρ₀: Background (hot gas) density
- ρ_contrast: Density contrast (cold gas density = ρ₀ * ρ_contrast)
- R: Cylinder radius
- δ: Smoothing thickness
- r: Distance from cylinder axis

### Velocity Profile
The shear velocity follows a similar tanh profile:
```
vₓ(r) = v_shear * (-tanh((r - R)/δ_vel))
```
Where:
- v_shear: Maximum shear velocity
- δ_vel: Velocity smoothing thickness

### Perturbations
Velocity perturbations are added in the transition region to seed the instability:
- **Magnitude**: v_pert * exp(-((r-R)/δ_vel)²) * f(x)
- **Direction**: Cylindrical coordinates (perpendicular to cylinder axis)
- **Spatial variation**: 
  - Sinusoidal: f(x) = sin(2π x/λ) when λ > 0
  - Localized: f(x) = exp(-((x-x₀)/δ)²) when λ = 0

## Parameters

### Required Parameters
- `rho_0`: Background density (hot gas)
- `pgas_0`: Background pressure (constant throughout domain)
- `density_contrast`: Density ratio (cold/hot)
- `radius`: Cylinder radius
- `vel_shear`: Shear velocity magnitude
- `smoothing_thickness`: Density transition thickness
- `gamma`: Adiabatic index

### Optional Parameters
- `smoothing_thickness_vel`: Velocity smoothing (default: same as density smoothing)
- `vel_pert`: Perturbation amplitude (default: 0.1)
- `lambda_pert`: Perturbation wavelength (default: 1.0, set to 0 for localized)
- `pert_loc`: Localized perturbation location (default: 0.0)
- `uniform_bfield`: Uniform magnetic field [Bx, By, Bz] (for MHD)
- `noisy_ic`: Add random noise (default: false)
- `random_seed`: Random seed (default: 12345)

## Comparison with Athena++ Reference

### Key Similarities
1. **Coordinate system**: Same cylindrical geometry with x as cylinder axis
2. **Density profile**: Identical tanh-based density distribution
3. **Velocity structure**: Same shear profile and perturbation approach
4. **Energy initialization**: Consistent total energy calculation including kinetic and thermal components

### Key Differences
1. **Code structure**: Enzo-E uses different field access patterns and class hierarchy
2. **Parameter system**: Enzo-E uses ParameterGroup instead of ParameterInput
3. **Field naming**: Enzo-E field names differ from Athena++ (e.g., "total_energy" vs "energy")
4. **Magnetic field handling**: Enzo-E has different MHD field structure

### Physics Equivalence
The physical setup is mathematically equivalent to the Athena++ implementation:
- Same density contrast and smoothing
- Identical velocity shear and perturbation profiles
- Consistent energy initialization
- Compatible magnetic field handling (when enabled)

## Usage Example

```cpp
Initial {
    list = ["kelvin_helmholtz"];
    
    kelvin_helmholtz {
        rho_0 = 1.0;                    # Background density
        pgas_0 = 1.0;                   # Background pressure
        density_contrast = 10.0;        # Cold/hot density ratio
        radius = 0.3;                   # Cylinder radius
        vel_shear = 1.0;                # Shear velocity
        vel_pert = 0.1;                 # Perturbation amplitude
        lambda_pert = 1.0;              # Perturbation wavelength
        smoothing_thickness = 0.05;     # Transition thickness
        gamma = 1.66667;                # Adiabatic index
    }
}
```

## Expected Behavior

1. **Initial state**: Cylindrical density and velocity structure
2. **Early evolution**: Growth of KH modes at the interface
3. **Nonlinear phase**: Formation of characteristic KH vortices
4. **Late evolution**: Turbulent mixing between cold and hot phases

The instability growth rate and characteristic scales depend on the density contrast, shear velocity, and smoothing thickness parameters.

## Files Modified/Created

1. **EnzoInitialKelvinHelmholtz.hpp**: Header file with class declaration
2. **EnzoInitialKelvinHelmholtz.cpp**: Implementation file
3. **initial.hpp**: Added include for new header
4. **CMakeLists.txt**: Added source files to build
5. **EnzoProblem.cpp**: Added factory method registration
6. **kelvin_helmholtz_example.in**: Example input file

This implementation provides a faithful reproduction of the Athena++ KH instability test case while being properly integrated into the Enzo-E framework.