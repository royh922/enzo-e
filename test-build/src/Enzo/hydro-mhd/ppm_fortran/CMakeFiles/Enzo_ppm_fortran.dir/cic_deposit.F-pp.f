# 1 "/Users/juitenghsu/enzo-e/src/Enzo/hydro-mhd/ppm_fortran/cic_deposit.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/Users/juitenghsu/enzo-e/src/Enzo/hydro-mhd/ppm_fortran/cic_deposit.F"


# 1 "/Users/juitenghsu/enzo-e/src/Enzo/fortran.h" 1

# 1 "/Users/juitenghsu/enzo-e/src/Enzo/enzo_defines.hpp" 1



























# 47 "/Users/juitenghsu/enzo-e/src/Enzo/enzo_defines.hpp"

# 2 "/Users/juitenghsu/enzo-e/src/Enzo/fortran.h" 2
















# 3 "/Users/juitenghsu/enzo-e/src/Enzo/hydro-mhd/ppm_fortran/cic_deposit.F" 2
c=======================================================================
c//////////////////////  SUBROUTINE CIC_DEPOSIT  \\\\\\\\\\\\\\\\\\\\\\c

      subroutine cic_deposit(posx, posy, posz, ndim, npositions, 
     &                      mass, field, leftedge, 
     &                      dim1, dim2, dim3, cellsize, cloudsize)

c
c  PERFORMS 1/2/3D CLOUD-IN-CELL INTERPOLATION FROM FIELD TO SUMFIELD
c
c  written by: Greg Bryan
c  date:       January, 1998
c  modified1:
c
c  PURPOSE: This routine performs a three-dimension, second-order
c           interpolation from field to sumfield (without clearing sumfield
c           first) at the positions specified.
c           The CIC cloud size is set by cloudsize, which should be
c           equal to or less than cellsize.  If cloudsize is equal to
c           cellsize, then this is the typical CIC deposit, but if it
c           is less than cellsize, then the effective cloud size is
c           smaller and the overlap calculations are adjusted.
c
c  INPUTS:
c     ndim       - dimensionality
c     cellsize   - the cell size of field
c     cloudsize  - size of the particle "cloud" (must be <= cellsize)
c     dim1,2,3   - real dimensions of field
c     leftedge   - the left edge(s) of field
c     npositions - number of particles
c     posx,y,z   - particle positions
c     sumfield   - 1D field (length npositions) of masses
c
c  OUTPUT ARGUMENTS: 
c     field      - field to be deposited to
c
c  EXTERNALS: 
c
c  LOCALS:
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

# 48 "/Users/juitenghsu/enzo-e/src/Enzo/hydro-mhd/ppm_fortran/cic_deposit.F" 2
c
c-----------------------------------------------------------------------
c
c  argument declarations
c
      integer (kind=4) dim1, dim2, dim3, npositions, ndim
      real(kind=8) posx(npositions), posy(npositions), posz(npositions),
     &        leftedge(3)
      real(kind=8)    mass(npositions), field(dim1, dim2, dim3), cellsize,
     &          cloudsize
c
c  locals
c
      integer (kind=4) iii, jjj, kkk
      integer (kind=4) i1, j1, k1, n
      real(kind=8)    xpos, ypos, zpos, dx, dy, dz
      real(kind=8) edge1, edge2, edge3, fact, shift, half, refine
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////
c=======================================================================
c

!     write(0,*) npositions, leftedge, dim1, dim2, dim3, cellsize

c      print*, 'cic_deposit.F mass(1)=',mass(1)
      fact = 1._PKIND/cellsize
      refine = cellsize/cloudsize
      half = 0.5001_PKIND/refine
      shift = 0.5_PKIND/refine
      edge1 = REAL(dim1,PKIND) - half
      edge2 = REAL(dim2,PKIND) - half
      edge3 = REAL(dim3,PKIND) - half
c
      if (cloudsize .gt. cellsize) then
         print*, 'ERROR'
C         ERROR_MESSAGE
         write(6,*) "cloudsize > cellsize in cic_deposit!"
         write(6,*) cloudsize, cellsize
      endif
c
c     1D
c
      if (ndim .eq. 1) then
c
         do n=1, npositions
c
c           Compute the position of the central cell
c
            xpos = min(max((posx(n) - leftedge(1))*fact, half), edge1)
c
c           Convert this into an integer (kind=4) index
c
            i1  = int(xpos - shift,IKIND) + 1
c
c           Compute the weights
c
            dx = min((REAL(i1,RKIND)+shift-xpos)*refine, 1.0_RKIND)
c
c           Interpolate from field into sumfield
c
            field(i1  ,1,1) = field(i1  ,1,1) + mass(n)*dx
            field(i1+1,1,1) = field(i1+1,1,1) + mass(n)*(1._RKIND-dx)
c
         enddo
c
      endif
c
c     2D
c
      if (ndim .eq. 2) then
c
         do n=1, npositions
c
c           Compute the position of the central cell
c
            xpos = min(max((posx(n) - leftedge(1))*fact, half), edge1)
            ypos = min(max((posy(n) - leftedge(2))*fact, half), edge2)
c
c           Convert this into an integer (kind=4) index
c
            i1  = int(xpos - shift,IKIND) + 1
            j1  = int(ypos - shift,IKIND) + 1
c
c           Compute the weights
c
            dx = min((REAL(i1,RKIND)+shift-xpos)*refine, 1.0_RKIND)
            dy = min((REAL(j1,RKIND)+shift-ypos)*refine, 1.0_RKIND)
c
c           Interpolate from field into sumfield
c
            field(i1  ,j1  ,1) = field(i1  ,j1  ,1) +
     &                           mass(n)*      dx *      dy
            field(i1+1,j1  ,1) = field(i1+1,j1  ,1) +
     &                           mass(n)*(1._RKIND-dx)*      dy
            field(i1  ,j1+1,1) = field(i1  ,j1+1,1) + 
     &                           mass(n)*      dx *(1._RKIND-dy)
            field(i1+1,j1+1,1) = field(i1+1,j1+1,1) +
     &                           mass(n)*(1._RKIND-dx)*(1._RKIND-dy)
c
         enddo
c
      endif
c
c     3D
c
      if (ndim .eq. 3) then
c     
         do n=1, npositions
c     
c     Compute the position of the central cell
c    
            xpos = min(max((posx(n) - leftedge(1))*fact, half),edge1)
            ypos = min(max((posy(n) - leftedge(2))*fact, half),edge2)
            zpos = min(max((posz(n) - leftedge(3))*fact, half),edge3)
c     
c     Convert this into an integer (kind=4) index
c     
            i1  = int(xpos - shift,IKIND) + 1
            j1  = int(ypos - shift,IKIND) + 1
            k1  = int(zpos - shift,IKIND) + 1

c     
c     Compute the weights
c     
            dx = min((REAL(i1,RKIND)+shift-xpos)*refine, 1.0_RKIND)
            dy = min((REAL(j1,RKIND)+shift-ypos)*refine, 1.0_RKIND)
            dz = min((REAL(k1,RKIND)+shift-zpos)*refine, 1.0_RKIND)
c     
c     Interpolate from field into sumfield
c
c            print*, 'DEBUG_CIC ',dx,dy,dz,i1,j1,k1,field(i1,j1,k1)
            field(i1  ,j1  ,k1  ) = field(i1  ,j1  ,k1  ) +
     &              mass(n)*     dx *     dy *    dz
            field(i1+1,j1  ,k1  ) = field(i1+1,j1  ,k1  ) +
     &              mass(n)*(1._RKIND-dx)*     dy *    dz
            field(i1  ,j1+1,k1  ) = field(i1  ,j1+1,k1  ) + 
     &              mass(n)*     dx *(1._RKIND-dy)*    dz
            field(i1+1,j1+1,k1  ) = field(i1+1,j1+1,k1  ) +
     &              mass(n)*(1._RKIND-dx)*(1._RKIND-dy)*    dz
            field(i1  ,j1  ,k1+1) = field(i1  ,j1  ,k1+1) +
     &              mass(n)*     dx *     dy *(1._RKIND-dz)
            field(i1+1,j1  ,k1+1) = field(i1+1,j1  ,k1+1) +
     &              mass(n)*(1._RKIND-dx)*     dy *(1._RKIND-dz)
            field(i1  ,j1+1,k1+1) = field(i1  ,j1+1,k1+1) + 
     &              mass(n)*     dx *(1._RKIND-dy)*(1._RKIND-dz)
            field(i1+1,j1+1,k1+1) = field(i1+1,j1+1,k1+1) +
     &           mass(n)*(1._RKIND-dx)*(1._RKIND-dy)*(1._RKIND-dz)

c     
         enddo
      endif
c
      return

      end
