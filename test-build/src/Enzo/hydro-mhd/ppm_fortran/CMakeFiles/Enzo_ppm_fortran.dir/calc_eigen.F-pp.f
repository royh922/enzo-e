# 1 "/Users/juitenghsu/enzo-e/src/Enzo/hydro-mhd/ppm_fortran/calc_eigen.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/Users/juitenghsu/enzo-e/src/Enzo/hydro-mhd/ppm_fortran/calc_eigen.F"

# 1 "/Users/juitenghsu/enzo-e/src/Enzo/fortran.h" 1

# 1 "/Users/juitenghsu/enzo-e/src/Enzo/enzo_defines.hpp" 1



























# 47 "/Users/juitenghsu/enzo-e/src/Enzo/enzo_defines.hpp"

# 2 "/Users/juitenghsu/enzo-e/src/Enzo/fortran.h" 2
















# 2 "/Users/juitenghsu/enzo-e/src/Enzo/hydro-mhd/ppm_fortran/calc_eigen.F" 2

c=======================================================================
c///////////////////////  SUBROUTINE CALC_EIGEN  \\\\\\\\\\\\\\\\\\\\\\\\c

      subroutine calc_eigen(dslice, cs, is, ie, idim, lem, rem)
c
c  COMPUTES THE LEFT AND RIGHT EIGENMATRICES FOR ADIABATIC HYDRODYNAMICS
c
c  written by: John Wise
c  date:       June, 2010
c  modified1:  
c
c  PURPOSE:  To project the primitive variables into the characteristic
c    variables and vice-versa, these matrices are needed.  We use these
c    variables to calculate the differences during interpolation.  This
c    makes the interpolation total variation diminishing (LeVeque 2002).
c
c  INPUT:
c    dslice - extracted 2d slice of the density, d
c    cs     - sound speed
c    
c  OUTPUT:
c    lem,rem  - left and right eigenmatrices
c
c-----------------------------------------------------------------------
c
      implicit NONE


# 1 "/Users/juitenghsu/enzo-e/src/Enzo/fortran_types.h" 1

# 20 "/Users/juitenghsu/enzo-e/src/Enzo/fortran_types.h"
      integer, parameter :: PKIND=8
      integer, parameter :: RKIND=8
# 33 "/Users/juitenghsu/enzo-e/src/Enzo/fortran_types.h"
      integer, parameter :: IKIND=4
      integer, parameter :: LKIND=4
# 74 "/Users/juitenghsu/enzo-e/src/Enzo/fortran_types.h"

# 31 "/Users/juitenghsu/enzo-e/src/Enzo/hydro-mhd/ppm_fortran/calc_eigen.F" 2
c
c     Arguments
c
      integer (kind=4) idim, is, ie
      real(kind=8) dslice(idim), cs(idim)
      real(kind=8) lem(idim,5,5), rem(idim,5,5)
c
c     Locals
c
      integer (kind=4) i, n, m
      real(kind=8) csi, csqi, csq
c
      do m = 1,5
         do n = 1,5
            do i = 1, idim
               lem(i,n,m) = 0._RKIND
               rem(i,n,m) = 0._RKIND
            enddo
         enddo
      enddo

      do i = 1, idim
         csi = 1._RKIND/cs(i)
         csq = cs(i)**2
         csqi = 1._RKIND/csq
c
c     Left eigenvectors
c         
         lem(i,2,1) = -0.5_RKIND*dslice(i)*csi
         lem(i,5,1) = 0.5_RKIND*csqi
         lem(i,1,2) = 1._RKIND
         lem(i,5,2) = -csqi
         lem(i,3,3) = 1._RKIND
         lem(i,4,4) = 1._RKIND
         lem(i,2,5) = -lem(i,2,1)
         lem(i,5,5) = lem(i,5,1)
c
c     Right eigenvectors
c
         rem(i,1,1) = 1._RKIND
         rem(i,1,2) = -cs(i)/dslice(i)
         rem(i,1,5) = csq
         rem(i,2,1) = 1._RKIND
         rem(i,3,3) = 1._RKIND
         rem(i,4,4) = 1._RKIND
         rem(i,5,1) = 1._RKIND
         rem(i,5,2) = -rem(i,1,2)
         rem(i,5,5) = csq

      enddo

      end
