# 1 "/Users/juitenghsu/enzo-e/src/Enzo/hydro-mhd/ppm_fortran/flux_twoshock.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/Users/juitenghsu/enzo-e/src/Enzo/hydro-mhd/ppm_fortran/flux_twoshock.F"


# 1 "/Users/juitenghsu/enzo-e/src/Enzo/fortran.h" 1

# 1 "/Users/juitenghsu/enzo-e/src/Enzo/enzo_defines.hpp" 1



























# 47 "/Users/juitenghsu/enzo-e/src/Enzo/enzo_defines.hpp"

# 2 "/Users/juitenghsu/enzo-e/src/Enzo/fortran.h" 2
















# 3 "/Users/juitenghsu/enzo-e/src/Enzo/hydro-mhd/ppm_fortran/flux_twoshock.F" 2
c=======================================================================
c//////////////////////  SUBROUTINE FLUX_TWOSHOCK \\\\\\\\\\\\\\\\\\\\\c

      subroutine flux_twoshock(
     $     dslice, eslice, geslice, uslice, vslice, wslice, dx,
     $     diffcoef, idim, jdim, i1, i2, j1, j2, dt, gamma, idiff,
     $     idual, eta1, ifallback,
     $     dls, drs, pls, prs, gels, gers, 
     $     uls, urs, vls, vrs, wls, wrs, pbar, ubar,
     $     df, ef, uf, vf, wf, gef, ges,
     $     ncolor, colslice, colls, colrs, colf, ierror
     $     )
c
c-----------------------------------------------------------------------
c
c  COMPUTE EULERIAN FLUXES FOR THE (LAGRANGIAN) TWO-SHOCK RIEMANN SOLVER
c
c  written by: Jim Stone
c  date:       January, 1991
c  modified1:  June, 1993 by Greg Bryan (Lagrange+remap --> Eulerian)
c  modified2:  July, 1994 by GB; switched to slicewise and put all the
c                                information into the argument line
c  modified3:  July, 1994 by GB; moved diffusion coefficient to own routine
c  modified4:  Febuary, 2005 by Alexei Kritsuk; fixed a bug in RAREFACTION1
c                            and a bug in update eq. 3.1 with idiff != 0.
c  modified5:  Sept. 2005 by DC, fixed the flux array to work with cosmology and AMR
c  modified6:  June, 2010 by John Wise; transferred from euler.src into its
c                         own routine
c
      implicit none

# 1 "/Users/juitenghsu/enzo-e/src/Enzo/fortran_types.h" 1

# 20 "/Users/juitenghsu/enzo-e/src/Enzo/fortran_types.h"
      integer, parameter :: PKIND=8
      integer, parameter :: RKIND=8
# 33 "/Users/juitenghsu/enzo-e/src/Enzo/fortran_types.h"
      integer, parameter :: IKIND=4
      integer, parameter :: LKIND=4
# 74 "/Users/juitenghsu/enzo-e/src/Enzo/fortran_types.h"

# 34 "/Users/juitenghsu/enzo-e/src/Enzo/hydro-mhd/ppm_fortran/flux_twoshock.F" 2
      integer (kind=4) ijkn, ierror
      parameter (ijkn=(1024+6))
c
c-----------------------------------------------------------------------
c     Arguments
c
      integer (kind=4) i1, i2, j1, j2, idim, jdim, ncolor, idiff, idual,
     $     ifallback
      real(kind=8) gamma, dt, eta1
      real(kind=8) ubar(idim,jdim), pbar(idim,jdim), dx(idim),
     $     uls(idim,jdim), vls(idim,jdim), wls(idim,jdim),
     $     urs(idim,jdim), vrs(idim,jdim), wrs(idim,jdim),
     $     pls(idim,jdim), dls(idim,jdim), gels(idim,jdim),
     $     prs(idim,jdim), drs(idim,jdim), gers(idim,jdim),
     $     colls(idim,jdim,ncolor), colrs(idim,jdim,ncolor),
     $     diffcoef(idim,jdim), dslice(idim,jdim), 
     $     uslice(idim,jdim), vslice(idim,jdim), wslice(idim,jdim),
     $     eslice(idim,jdim), geslice(idim,jdim),
     $     colslice(idim,jdim,ncolor)
      real(kind=8) df(idim,jdim), ef(idim,jdim), gef(idim,jdim),
     $     uf(idim,jdim), vf(idim,jdim), wf(idim,jdim),
     $     colf(idim,jdim,ncolor), ges(idim,jdim)

c
c-----------------------------------------------------------------------
c     Locals
c
      integer (kind=4) i, j, n
      real(kind=8) qa, qb, qc, frac, one
      real(kind=8) sn(ijkn), u0(ijkn), p0(ijkn), d0(ijkn), c0(ijkn), z0(ijkn),
     $     dbar(ijkn), cbar(ijkn), l0(ijkn), lbar(ijkn), 
     $     ub(ijkn), vb(ijkn), wb(ijkn), db(ijkn), pb(ijkn), 
     $     cb(ijkn), eb(ijkn), geb(ijkn), colb(ijkn,20),
     $     upb(ijkn), dub(ijkn), duub(ijkn), duvb(ijkn), duwb(ijkn),
     $     dueb(ijkn), dugeb(ijkn), pcent(ijkn)
      parameter (one=1._RKIND)
c
c-----------------------------------------------------------------------
c
      qa = (gamma + 1._RKIND)/(2._RKIND*gamma)
      qb = (gamma - 1._RKIND)/(gamma + 1._RKIND)
c
c  Loop over sweep lines (in this slice)
c
      do j=j1, j2
c
c  Evaluate time-averaged quantities 
c   (see Colella, Siam J Sci Stat Comput 1982, 3, 77.  Appendix)
c
       do i=i1, i2+1
c
          sn(i)    = sign(one, -ubar(i,j))
c
c  Collect values of interest depending on which way fluid is flowing
c     
          if (sn(i) .lt. 0._RKIND) then
             u0(i) = uls(i,j)
             p0(i) = pls(i,j)
             d0(i) = dls(i,j)
c             vb(i) = vls(i,j)
c             wb(i) = wls(i,j)
c             geb(i) = gels(i,j)
          else
             u0(i) = urs(i,j)
             p0(i) = prs(i,j)
             d0(i) = drs(i,j)
c             vb(i) = vrs(i,j)
c             wb(i) = wrs(i,j)
c             geb(i) = gers(i,j)
          endif
          c0(i)    = sqrt(max(gamma*p0(i)/d0(i), 1.d-20))
          z0(i)    = c0(i)*d0(i)*sqrt(max(
     &               1._RKIND + qa*(pbar(i,j)/p0(i)-1._RKIND), 1.d-20))
c          write(24,*) i, u0(i), uls(i,j), urs(i,j), c0(i)
       enddo
c
c      Repeat for color variables (now moved below)
c
# 123 "/Users/juitenghsu/enzo-e/src/Enzo/hydro-mhd/ppm_fortran/flux_twoshock.F"
c
c  Compute equivalent bar (inside shock & rarefaction) values for density
c    and sound speed
c
       do i=i1, i2+1
          dbar(i)  = 1._RKIND/(1._RKIND/d0(i) - (pbar(i,j)-p0(i))/
     &         max(z0(i)**2, 1.d-20))
          cbar(i)  = sqrt(max(gamma*pbar(i,j)/dbar(i), 1.d-20))
       enddo
c
c
c  Find lambda values for the shock and rarefaction
c
       do i=i1, i2+1
         if (pbar(i,j) .lt. p0(i)) then
            l0(i)   = u0(i)*sn(i) + c0(i)
            lbar(i) = sn(i)*ubar(i,j)+cbar(i)
         else
            l0(i)   = u0(i)*sn(i) + z0(i)/d0(i)
c            l0(i)   = ubar(i,j)*sn(i) + z0(i)/d0(i)
            lbar(i) = l0(i)
         endif
c
c  Compute values for inside a rarefaction fan
c     (method described in Colella, 1982)
c

c
# 162 "/Users/juitenghsu/enzo-e/src/Enzo/hydro-mhd/ppm_fortran/flux_twoshock.F"
c
# 177 "/Users/juitenghsu/enzo-e/src/Enzo/hydro-mhd/ppm_fortran/flux_twoshock.F"
c

c
c  Compute values for inside a rarefaction fan
c     (linear interpolation between end states, as suggested in PPM ref)
c
         frac     = l0(i) - lbar(i)
         if (frac .lt. 1.d-20) frac = 1.d-20
         frac     = (0._RKIND - lbar(i))/frac
         frac     = min(max(frac, 0._RKIND), 1._RKIND)
         pb(i)    = p0(i)*frac + pbar(i,j)*(1._RKIND - frac)
         db(i)    = d0(i)*frac + dbar(i  )*(1._RKIND - frac)
         ub(i)    = u0(i)*frac + ubar(i,j)*(1._RKIND - frac)
c

c
       enddo
c
c  Cull appropriate states depending on where eulerian position is in solution
c    (lbar >= 0 --> inside post-shock region,
c     l0   <  0 --> outside shock/rarefaction wave,
c     otherwise --> inside rarefaction wave).
c
       do i=i1, i2+1
         if (lbar(i) .ge. 0._RKIND) then !AK
            pb(i) = pbar(i,j)
            db(i) = dbar(i  )
            ub(i) = ubar(i,j)
         endif
         if (l0  (i) .lt. 0._RKIND) then !AK
            pb(i) = p0(i)
            db(i) = d0(i)
            ub(i) = u0(i)
         endif
c         write(22,*) i, ub(i), u0(i), ubar(i,j), sn(i)
c
c         if (db(i) .gt. dbar(i) .and. db(i) .gt. d0(i))
c     &      write(6,1000) i, j, db(i), dbar(i), d0(i)
c 1000  format('euler: ', 2(i4), 3(e12.4))
c
       enddo
c
# 233 "/Users/juitenghsu/enzo-e/src/Enzo/hydro-mhd/ppm_fortran/flux_twoshock.F"
c
c  Collect values of interest depending on which way fluid is flowing
c     (note: new placement for this - uses ub instead of ubar for consistency)
c
       do i=i1, i2+1
          if (ub(i) .gt. 0._RKIND) then !AK
             vb(i) = vls(i,j)
             wb(i) = wls(i,j)
             geb(i) = gels(i,j)
          else
             vb(i) = vrs(i,j)
             wb(i) = wrs(i,j)
             geb(i) = gers(i,j)
          endif
       enddo
c
c      Repeat for color variables
c
       do n=1,ncolor
          do i=i1, i2+1
             if (ub(i) .gt. 0._RKIND) then
                colb(i,n) = colls(i,j,n) * db(i)/dls(i,j)
             else
                colb(i,n) = colrs(i,j,n) * db(i)/drs(i,j)
             endif
          enddo
       enddo
c
c  Dual energy formalism: if sound speed squared is less than eta1*v^2 
c    then discard pbar,dbar,ubar in favour of p0,d0,u0.  This amounts
c     to assuming that we are outside the shocked region but the flow is
c     hypersonic so this should be true.  This is inserted because the
c     mechanism above occasionally fails in hypersonic flows.
c
# 279 "/Users/juitenghsu/enzo-e/src/Enzo/hydro-mhd/ppm_fortran/flux_twoshock.F"
c
c  Calculate total specific energy corresponding to this state
c     (and the specific gas energy).
c
       do i=i1, i2+1
         eb(i) = pb(i)/((gamma-1._RKIND)*db(i)) + 
     &           0.5_RKIND*(ub(i)**2 + vb(i)**2 + wb(i)**2)
c         geb(i) = pb(i)/((gamma-1._RKIND)*db(i))
       enddo
c
c  Compute terms in differenced hydro equations (eq. 3.1)
c
       if (idiff .ne. 0) then
c
c     ...with diffusion
c
         do i=i1,i2+1
            upb(i)  =  pb(i)*ub(i)
c
            dub(i)  =  ub(i)*db(i)  !AK
c
            duub(i) = dub(i)*ub(i) + diffcoef(i,j)*
     &           (dslice(i-1,j)*uslice(i-1,j) - dslice(i,j)*uslice(i,j))
c
            duvb(i) = dub(i)*vb(i)
            duwb(i) = dub(i)*wb(i)
c
c     (should we add diffusion to the cross velocities?  I doubt it)
c     I don't. This doubt kills Noh test problem in 2D at high resolution.
c     Diffusion has to be added to cross velocities. !AK        May 2005.
c
            duvb(i) = dub(i)*vb(i) + diffcoef(i,j)*
     &           (dslice(i-1,j)*vslice(i-1,j) - dslice(i,j)*vslice(i,j))
            duwb(i) = dub(i)*wb(i) + diffcoef(i,j)*
     &           (dslice(i-1,j)*wslice(i-1,j) - dslice(i,j)*wslice(i,j))
c
            dueb(i) = dub(i)*eb(i) + diffcoef(i,j)*
     &           (dslice(i-1,j)*eslice(i-1,j) - dslice(i,j)*eslice(i,j))
c
c     This update must be the last !AK
c
            dub(i)  =  dub(i) + diffcoef(i,j)*
     &           (dslice(i-1,j)               - dslice(i,j)            )
         enddo
c
c        If using dual energy formalism, compute dugeb
c
         if (idual .eq. 1) then
            do i=i1,i2+1
              dugeb(i) = dub(i)*geb(i) + diffcoef(i,j)*
     &         (dslice(i-1,j)*geslice(i-1,j) - dslice(i,j)*geslice(i,j))
            enddo
         endif
c
       else
c
c     ...and without
c
         do i=i1, i2+1
            upb(i)  =  pb(i)*ub(i)
            dub(i)  =  ub(i)*db(i)
            duub(i) = dub(i)*ub(i)
            duvb(i) = dub(i)*vb(i)
            duwb(i) = dub(i)*wb(i)
            dueb(i) = dub(i)*eb(i)
         enddo
c
         if (idual .eq. 1) then
            do i=i1, i2+1
               dugeb(i) = dub(i)*geb(i)
            enddo
         endif
c
       endif
c
c  Copy into flux slices (to return to caller)
c

!      do i=i1, i2+1
!         df(i,j) = dt*dub(i)
!         ef(i,j) = dt*(dueb(i) + upb(i))
!         uf(i,j) = dt*(duub(i) + pb(i))
!         vf(i,j) = dt*duvb(i)
!         wf(i,j) = dt*duwb(i)
!      enddo

       do i=i1, i2+1
          qc = dt/dx(i)
          df(i,j) = qc*dub(i)
          ef(i,j) = qc*(dueb(i) + upb(i))
          uf(i,j) = qc*(duub(i) + pb(i))
          vf(i,j) = qc*duvb(i)
          wf(i,j) = qc*duwb(i)
       enddo

c
       do n=1,ncolor
          do i=i1, i2+1
c             colf(i,j,n) = dt*dub(i)*colb(i,n)   ! color*dens conserved
             colf(i,j,n) = dt*ub(i)*colb(i,n)     ! color      conserved
          enddo
       enddo
c
c      Do the same for the gas energy if using the dual energy formalism
c         (note that we do not include the source term)
c      JHW (Jul 2010): Moved the source term to here.
c
       if (idual .eq. 1) then
          do i=i1, i2+1
             pcent(i) = max((gamma-1._RKIND)*geslice(i,j)*dslice(i,j), 
     $            1.d-20)
          enddo
          do i=i1, i2+1
             qc = dt/dx(i)
             gef(i,j) = qc*dugeb(i)
             ges(i,j) = qc * pcent(i) * (ub(i) - ub(i+1))
          enddo
       endif

c
c     Check here for negative densities and energies
c
       do i=i1, i2
          if (dslice(i,j) + (df(i,j)-df(i+1,j)) .le. 0._RKIND .or.
     $         eslice(i,j) .lt. 0._RKIND) then
             if (eslice(i,j) .lt. 0._RKIND) then
                write (6,*) 'flux_twoshock: eslice < 0: ', i, j
                write (6,*) qc*(dueb(i)-dueb(i+1)),
     &                     qc*(upb(i)-upb(i+1)), eslice(i,j)*dslice(i,j)
                write (6,*) eb(i-1), eb(i), eb(i+1)
                write (6,*) dueb(i-1), dueb(i), dueb(i+1)
                write (6,*) upb(i-1), upb(i), upb(i+1)
                ierror = 10303
             else
                write (6,*) 'flux_twoshock: dnu <= 0:', i, j
                ierror = 10302
             endif
             write (6,*) 'a', df(i,j)-df(i+1,j), dslice(i,j), 
     $            dt/dx(i), dt, lbar(i),l0(i)
             write (6,*) 'b', db(i-1), db(i), db(i+1)
             write (6,*) 'c', ub(i-1), ub(i), ub(i+1)
             write (6,*) 'd', pb(i-1), pb(i), pb(i+1)
             write (6,*) 'e', d0(i-1), d0(i), d0(i+1)
             write (6,*) 'f', u0(i-1), u0(i), u0(i+1)
             write (6,*) 'g', p0(i-1), p0(i), p0(i+1)
             write (6,*) 'h', ubar(i-1,j), ubar(i,j), ubar(i+1,j)
             write (6,*) 'i', pbar(i-1,j), pbar(i,j), pbar(i+1,j)
             write (6,*) 'j', dls(i-1,j), dls(i,j), dls(i+1,j)
             write (6,*) 'k', drs(i-1,j), drs(i,j), drs(i+1,j)
             write (6,*) 'l', uls(i-1,j), uls(i,j), uls(i+1,j)
             write (6,*) 'm', urs(i-1,j), urs(i,j), urs(i+1,j)
             write (6,*) 'n', pls(i-1,j), pls(i,j), pls(i+1,j)
             write (6,*) 'o', prs(i-1,j), prs(i,j), prs(i+1,j)
             write (6,*) 'p', dslice(i-1,j), dslice(i,j), dslice(i+1,j)
             write (6,*) 'q', uslice(i-1,j), uslice(i,j), uslice(i+1,j)
             write (6,*) 'r', geslice(i-1,j), geslice(i,j), 
     $            geslice(i+1,j)
             write (6,*) 's', eslice(i-1,j), eslice(i,j), eslice(i+1,j)
             write (6,*) 't', geslice(i-1,j)*dslice(i-1,j)*(gamma-1.0), 
     &                   geslice(i  ,j)*dslice(i  ,j)*(gamma-1.0), 
     &                   geslice(i+1,j)*dslice(i+1,j)*(gamma-1.0)
c             if (gravity .eq. 1) 
c     &        write (6,*) grslice(i-1,j), grslice(i,j), grslice(i+1,j)
             if (ifallback.eq.0) then
                write(0,*) 'stop in euler with e < 0'
                ierror = 10304
                print*, 'ERROR_MESSAGE'
             else
                print*, 'WARNING'
                write(0,*) 'Falling back to a more diffusive Riemann',
     $               'solver, HLL'
                call flux_hll(dslice, eslice, geslice, 
     $               uslice, vslice, wslice,
     $               dx, diffcoef, idim, jdim, i, i, j, j, dt, gamma,
     $               idiff, idual, eta1, ifallback,
     $               dls, drs, pls, prs, uls, urs, vls, vrs, wls, wrs,
     $               gels, gers, df, uf, vf, wf, ef, gef, ges, 
     $               ncolor, colslice, colls, colrs, colf, ierror)
             endif
          endif
       enddo
c
c     Now check for negative colors (Warning only)
c
# 490 "/Users/juitenghsu/enzo-e/src/Enzo/hydro-mhd/ppm_fortran/flux_twoshock.F"

      enddo                     ! j

      end
