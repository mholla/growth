c...  ------------------------------------------------------------------
      subroutine sdvini(statev,coords,nstatv,ncrds,noel,npt,layer,kspt)
c...  ------------------------------------------------------------------
      include 'aba_param.inc'

      dimension statev(nstatv)

      statev(1)=1.0d0
      
      return
      end

c...  ------------------------------------------------------------------
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     #rpl,ddsddt,drplde,drpldt,
     #stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     #ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     #celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
c...  ------------------------------------------------------------------
      include 'aba_param.inc'

      character*80 cmname
      dimension stress(ntens),statev(nstatv),
     # ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     # stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     # props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      call umat_neohooke(stress,statev,ddsdde,sse,
     #                   time,dtime,coords,props,dfgrd1,
     #                   ntens,ndi,nshr,nstatv,nprops,
     #                   noel,npt,kstep,kinc)
      
      return
      end

c...  ------------------------------------------------------------------
      subroutine umat_neohooke(stress,statev,ddsdde,sse,
     #                         time,dtime,coords,props,dfgrd1,
     #                         ntens,ndi,nshr,nstatv,nprops,
     #                         noel,npt,kstep,kinc))
c...  ------------------------------------------------------------------
c * " This UMAT is written for a Neohookean constitutive material.
c * " It differs from the built-in Abaqus model and may not match.

c * " Written by Maria Holland in Aug 2012
c * " Updated by Maria Holland in Sep 2017
c...  ------------------------------------------------------------------

      implicit none

c...  variables to be defined
      real*8  stress(ntens), ddsdde(ntens,ntens), statev(nstatv), sse

c...  variables passed in for information
      real*8  time(2), dtime, coords(3), props(nprops), dfgrd1(3,3)
      integer ntens, ndi, nshr, nstatv, nprops, noel, npt, kstep, kinc

c...  material properties
      real*8  lam,mu

c...  local variables      
      integer i,j
      real*8  detf, lnJ
      real*8  b(6), xi(6)
      real*8  sse, statev
      
      data xi/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/

c...  initialize material parameters 
      lam = props(1)              ! lame constant
      mu  = props(2)              ! shear modulus      

c...  calculate determinant of deformation gradient
      detf=+dfgrd1(1,1) * (dfgrd1(2,2)*dfgrd1(3,3)-dfgrd1(2,3)*dfgrd1(3,2))
     #     -dfgrd1(1,2) * (dfgrd1(2,1)*dfgrd1(3,3)-dfgrd1(2,3)*dfgrd1(3,1))
     #     +dfgrd1(1,3) * (dfgrd1(2,1)*dfgrd1(3,2)-dfgrd1(2,2)*dfgrd1(3,1))
     
      lnJ = dlog(detf)
     
c...  calculate left cauchy-green deformation tensor b = f * f^t
      b(1) = dfgrd1(1,1)*dfgrd1(1,1) + dfgrd1(1,2)*dfgrd1(1,2) + dfgrd1(1,3)*dfgrd1(1,3)
      b(2) = dfgrd1(2,1)*dfgrd1(2,1) + dfgrd1(2,2)*dfgrd1(2,2) + dfgrd1(2,3)*dfgrd1(2,3)
      b(3) = dfgrd1(3,1)*dfgrd1(3,1) + dfgrd1(3,2)*dfgrd1(3,2) + dfgrd1(3,3)*dfgrd1(3,3)
      b(4) = dfgrd1(1,1)*dfgrd1(2,1) + dfgrd1(1,2)*dfgrd1(2,2) + dfgrd1(1,3)*dfgrd1(2,3)
      b(5) = dfgrd1(1,1)*dfgrd1(3,1) + dfgrd1(1,2)*dfgrd1(3,2) + dfgrd1(1,3)*dfgrd1(3,3)
      b(6) = dfgrd1(2,1)*dfgrd1(3,1) + dfgrd1(2,2)*dfgrd1(3,2) + dfgrd1(2,3)*dfgrd1(3,3)

c...  calculate Cauchy stress
      do i = 1,ntens
        stress(i) = ((lam*lnJ-mu)*xi(i) + mu*b(i))/detf
      enddo      
      
c...  calculateelastic and geometric tangent
      ddsdde(1,1) = (lam + 2.d0*(lam*lnJe - mu))/detf + 2.d0*stress(1)
      ddsdde(2,2) = (lam + 2.d0*(lam*lnJe - mu))/detf + 2.d0*stress(2)
      ddsdde(3,3) = (lam + 2.d0*(lam*lnJe - mu))/detf + 2.d0*stress(3)
      ddsdde(1,2) = (lam)/detf
      ddsdde(1,3) = (lam)/detf
      ddsdde(2,3) = (lam)/detf
      ddsdde(1,4) = stress(4)
      ddsdde(2,4) = stress(4)
      ddsdde(3,4) = 0.d0
      ddsdde(4,4) = (lam*lnJe - mu)/detf + (stress(1) + stress(2))/2.d0
      
      if (ntens.eq.6) then
        ddsdde(1,5) = stress(5)
        ddsdde(2,5) = 0.d0
        ddsdde(3,5) = stress(5)
        ddsdde(1,6) = 0.d0
        ddsdde(2,6) = stress(6)
        ddsdde(3,6) = stress(6)
        ddsdde(5,5) = (lam*lnJe - mu)/detf + (stress(1) + stress(3))/2.d0
        ddsdde(6,6) = (lam*lnJe - mu)/detf + (stress(2) + stress(3))/2.d0
        ddsdde(4,5) = stress(6)/2.d0
        ddsdde(4,6) = stress(5)/2.d0
        ddsdde(5,6) = stress(4)/2.d0
      endif
    
c...  use symmetry to fill in the rest
      do i = 2,ntens
        do j = 1,i-1
          ddsdde(i,j) = ddsdde(j,i)
        enddo
      enddo    
      
c...  calculate strain energy
      sse = (lam*lnJ**2.d0 + mu*(b(1)+b(2)+b(3) - 3.d0 - 2.d0*lnJ))/2.d0

      return
      end
