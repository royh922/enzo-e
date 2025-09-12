# 1 "/Users/juitenghsu/enzo-e/src/Enzo/cosmology/expand_terms.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/Users/juitenghsu/enzo-e/src/Enzo/cosmology/expand_terms.F"

# 1 "/Users/juitenghsu/enzo-e/src/Enzo/fortran.h" 1

# 1 "/Users/juitenghsu/enzo-e/src/Enzo/enzo_defines.hpp" 1



























# 47 "/Users/juitenghsu/enzo-e/src/Enzo/enzo_defines.hpp"

# 2 "/Users/juitenghsu/enzo-e/src/Enzo/fortran.h" 2
















# 2 "/Users/juitenghsu/enzo-e/src/Enzo/cosmology/expand_terms.F" 2
c=======================================================================
c///////////////////////  SUBROUTINE EXPAND_TERMS  \\\\\\\\\\\\\\\\\\\\c

      subroutine expand_terms(rank, isize, idual, coef, imethod, gamma,
     &                        p, d, e, ge, u, v, w,
     &                        dold, eold, geold, uold, vold, wold,
     &                        icr, ecr, ecrold)
c
c  ADDS THE COMOVING EXPANSION TERMS TO THE PHYSICAL VARIABLES
c
c     written by: Greg Bryan
c     date:       February, 1996
c     modified1:
c
c  PURPOSE:
c         Note: p is modified on exit
c
c  INPUTS:
c    isize   - size of fields (treat as 1d)
c    idual   - dual energy flag (1 = on, 0 = off)
c    coef    - coefficent (dt * adot / a)
c    d       - density field
c    p       - pressure field (from total energy - 0.5v^2)
c    e,ge    - total energy and gas energy (specific)
c    u,v,w   - velocities
c
c  OUTPUTS:
c    d,e,ge,u,v,w - output fields
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

# 38 "/Users/juitenghsu/enzo-e/src/Enzo/cosmology/expand_terms.F" 2
c
c     Arguments
c
      integer isize, idual, imethod, rank, icr
      real (kind=8) gamma, coef
      real (kind=8) d(isize), p(isize), e(isize), ge(isize),
     &        u(isize), v(isize), w(isize), ecr(isize)
      real (kind=8) dold(isize), eold(isize), geold(isize),
     &        uold(isize), vold(isize), wold(isize), ecrold(isize)
c
c     Locals
c
      integer i
c
c\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\///////////////////////////////////
c=======================================================================
c
c   METHOD1 is sortof time-centered
c   METHOD2 is time-backward
c   METHOD3 is semi-implicit
c
c   (If this is changed, you must also change the pressure computation
c    in Grid_ComovingExpandsionTerms.C)
c

c
c     Do gas energy first (if necessary); the term is 3*p/d
c
      if (idual .eq. 1) then
         do i = 1, isize
# 76 "/Users/juitenghsu/enzo-e/src/Enzo/cosmology/expand_terms.F"
            ge(i) = ge(i)*(1._RKIND-coef)/(1._RKIND+coef)
c            if (ge(i)-coef*(3._RKIND-2._RKIND/(gamma-1._RKIND))*p(i)/d(i).le.0)
c     &         write(6,*) 'get:',i,ge(i),p(i),d(i)
c
c  this line should be there if gamma != 5/3:
c            ge(i) = ge(i) - coef*(3._RKIND - 2._RKIND/(gamma-1._RKIND))*p(i)/d(i)
c

         enddo
      endif
c
c     Now do total energy; the term is 3*p/d + v^2
c       (for zeus method (imethod=2), only use 3*p/d term)
c
      if (imethod .eq. 2) then
         if (icr .gt. 0) then
c           For CR, use implicit on both
            do i = 1, isize
               e(i) = e(i) * (1.0 - coef)/(1.0 + coef)
               ecr(i) = ecr(i) * (2.0 - coef)/(2.0 + coef)
            enddo
         else
            do i = 1, isize
               e(i) = max(e(i) - coef*6._RKIND*p(i)/(d(i)+dold(i)), 
     &              0.5_RKIND*e(i))
            enddo
         endif
      else
         do i = 1, isize
# 124 "/Users/juitenghsu/enzo-e/src/Enzo/cosmology/expand_terms.F"
            e(i) = e(i)*(1._RKIND-coef)/(1._RKIND+coef)
            e(i) = max(e(i) - coef*(3._RKIND - 
     &           2._RKIND/(gamma-1._RKIND))*p(i)/d(i), 0.5_RKIND*e(i))

         enddo
      endif
c
c     Velocity terms
c

c
c        i) sortof time-centered: */
c
# 147 "/Users/juitenghsu/enzo-e/src/Enzo/cosmology/expand_terms.F"
c
c        iii) time-forward */
c







c
c        iii) semi-implicit way: */
c

      do i = 1, isize
         u(i) = u(i)*(1._RKIND-0.5_RKIND*coef) / 
     &        (1._RKIND + 0.5_RKIND*coef)
         if (rank .gt. 1) v(i) = v(i)*(1._RKIND-0.5_RKIND*coef) / 
     &        (1._RKIND + 0.5_RKIND*coef)
         if (rank .gt. 2) w(i) = w(i)*(1._RKIND-0.5_RKIND*coef) / 
     &        (1._RKIND + 0.5_RKIND*coef)
      enddo

c
c
      return
      end
