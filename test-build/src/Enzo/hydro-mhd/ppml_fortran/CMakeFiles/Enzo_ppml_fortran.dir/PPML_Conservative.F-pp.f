# 1 "/Users/juitenghsu/enzo-e/src/Enzo/hydro-mhd/ppml_fortran/PPML_Conservative.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/Users/juitenghsu/enzo-e/src/Enzo/hydro-mhd/ppml_fortran/PPML_Conservative.F"
c     See LICENSE_PPML file for license and copyright information


# 1 "/Users/juitenghsu/enzo-e/src/Enzo/fortran.h" 1

# 1 "/Users/juitenghsu/enzo-e/src/Enzo/enzo_defines.hpp" 1



























# 47 "/Users/juitenghsu/enzo-e/src/Enzo/enzo_defines.hpp"

# 2 "/Users/juitenghsu/enzo-e/src/Enzo/fortran.h" 2
















# 4 "/Users/juitenghsu/enzo-e/src/Enzo/hydro-mhd/ppml_fortran/PPML_Conservative.F" 2

      Subroutine Conserv(nx,ny,nz,dn,vx,vy,vz,bx,by,bz,
     &                   qu1,qu2,qu3,qu4,qu5,qu6,qu7)
         Implicit NONE

      Integer nx,ny,nz,i,j,k
      real (kind=8) dn(nx,ny,nz)           
      real (kind=8) vx(nx,ny,nz),vy(nx,ny,nz),vz(nx,ny,nz)           
      real (kind=8) bx(nx,ny,nz),by(nx,ny,nz),bz(nx,ny,nz)           

      real (kind=8) qu1(nx,ny,nz)           
      real (kind=8) qu2(nx,ny,nz),qu3(nx,ny,nz),qu4(nx,ny,nz)           
      real (kind=8) qu5(nx,ny,nz),qu6(nx,ny,nz),qu7(nx,ny,nz)           
          
      Do k=1,nz
       Do j=1,ny
        Do i=1,nx
             QU1(i,j,k)=dn(i,j,k)
             QU2(i,j,k)=dn(i,j,k)*vx(i,j,k)
             QU3(i,j,k)=dn(i,j,k)*vy(i,j,k)
             QU4(i,j,k)=dn(i,j,k)*vz(i,j,k)
             QU5(i,j,k)=bx(i,j,k)
             QU6(i,j,k)=by(i,j,k)
             QU7(i,j,k)=bz(i,j,k)
                Enddo
           Enddo        
      Enddo
           
         RETURN
         END
