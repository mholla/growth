c...  ------------------------------------------------------------------
      subroutine sdvini(statev,coords,nstatv,ncrds,noel,npt,layer,kspt)
c...  ------------------------------------------------------------------
      include 'aba_param.inc'

      dimension statev(nstatv)

      statev(1)=1.0d0    ! the_g
      statev(2)=1.0d0    ! J_g
      statev(3)=1.0d0    ! J_e      
      
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

      call umat_fiber_morph_Abaqus(stress,statev,ddsdde,sse,
     #                             time,dtime,coords,props,dfgrd1,
     #                             ntens,ndi,nshr,nstatv,nprops,
     #                             noel,npt,kstep,kinc)

      return
      end

c...  ------------------------------------------------------------------
      subroutine umat_fiber_morph_Abaqus(stress,statev,ddsdde,sse,
     #                                   time,dtime,coords,props,dfgrd1,
     #                                   ntens,ndi,nshr,nstatv,nprops,
     #                                   noel,npt,kstep,kinc)
c...  ------------------------------------------------------------------
c * " This UMAT is written for the Abaqus Neohookean constitutive material.
c * " It implements morphogenetic fiber growth in the direction of xn0.
c * " It is not suitable for use with local material directions.
c * " The growth rate can be easily changed (material properties and update formula)

c * " Written by Maria Holland in Aug 2012
c * " Updated by Maria Holland in Sep 2017
c...  ------------------------------------------------------------------
      
      implicit none

c...  variables to be defined
      real*8  stress(ntens),statev(nstatv),ddsdde(ntens,ntens),sse

c...  variables passed in for information
      real*8  time(2),dtime,coords(3),props(nprops),dfgrd1(3,3)
      integer ntens,ndi,nshr,nstatv,nprops,noel,npt,kstep,kinc

c...  material properties
      real*8  lam, mu, tmax, tau, xn0(3), k

c...  local variables      
      integer i,j
      real*8  norm, xn(3)
      real*8  theg, detf, detfe, I13e, J23e, vol
      real*8  fe(3,3), be(6), xi(6)
      
      data xi/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/

c...  initialize material parameters 
      lam    = props(1)              ! lame constant
      mu     = props(2)              ! shear modulus
      tmax   = props(3)
      tau    = props(4)
      xn0(1) = props(5)              ! f0[1]
      xn0(2) = props(6)              ! f0[2]
      xn0(3) = props(7)              ! f0[3]

      k = (3.d0*lam + 2.d0*mu)/3.d0
           
c...  calculate deformed elastic normal
      norm = sqrt(xn0(1)*xn0(1) + xn0(2)*xn0(2) + xn0(3)*xn0(3))
      xn(1) = (dfgrd1(1,1)*xn0(1) + dfgrd1(1,2)*xn0(2) + dfgrd1(1,3)*xn0(3))/norm
      xn(2) = (dfgrd1(2,1)*xn0(1) + dfgrd1(2,2)*xn0(2) + dfgrd1(2,3)*xn0(3))/norm
      xn(3) = (dfgrd1(3,1)*xn0(1) + dfgrd1(3,2)*xn0(2) + dfgrd1(3,3)*xn0(3))/norm

c...  update growth factor
      theg = (tmax-1.d0)*(1.d0-exp(-time(2)/tau)) + 1.d0
      statev(1) = theg
      statev(2) = theg     
      
c...  calculate elastic tensor Fe = F + ((1-theg)/theg)*(f \otimes f0)
      fe(1,1) = dfgrd1(1,1) + ((1.d0-theg)/theg)*(xn(1)*xn0(1))
      fe(1,2) = dfgrd1(1,2) + ((1.d0-theg)/theg)*(xn(1)*xn0(2))
      fe(1,3) = dfgrd1(1,3) + ((1.d0-theg)/theg)*(xn(1)*xn0(3))
      fe(2,1) = dfgrd1(2,1) + ((1.d0-theg)/theg)*(xn(2)*xn0(1))
      fe(2,2) = dfgrd1(2,2) + ((1.d0-theg)/theg)*(xn(2)*xn0(2)) 
      fe(2,3) = dfgrd1(2,3) + ((1.d0-theg)/theg)*(xn(2)*xn0(3))
      fe(3,1) = dfgrd1(3,1) + ((1.d0-theg)/theg)*(xn(3)*xn0(1))
      fe(3,2) = dfgrd1(3,2) + ((1.d0-theg)/theg)*(xn(3)*xn0(2))
      fe(3,3) = dfgrd1(3,3) + ((1.d0-theg)/theg)*(xn(3)*xn0(3))

c...  calculate determinant of deformation gradient
      detf=+dfgrd1(1,1) * (dfgrd1(2,2)*dfgrd1(3,3)-dfgrd1(2,3)*dfgrd1(3,2))
     #     -dfgrd1(1,2) * (dfgrd1(2,1)*dfgrd1(3,3)-dfgrd1(2,3)*dfgrd1(3,1))
     #     +dfgrd1(1,3) * (dfgrd1(2,1)*dfgrd1(3,2)-dfgrd1(2,2)*dfgrd1(3,1))
     
      detfe = +fe(1,1) * (fe(2,2)*fe(3,3)-fe(2,3)*fe(3,2))
     #        -fe(1,2) * (fe(2,1)*fe(3,3)-fe(2,3)*fe(3,1))
     #        +fe(1,3) * (fe(2,1)*fe(3,2)-fe(2,2)*fe(3,1))
     
      J23e = detfe**(-2.d0/3.d0)
      statev(3) = detfe
        
c...  calculate elastic left cauchy-green deformation tensor be = fe * fe^t
      be(1) = fe(1,1)*fe(1,1) + fe(1,2)*fe(1,2) + fe(1,3)*fe(1,3)
      be(2) = fe(2,1)*fe(2,1) + fe(2,2)*fe(2,2) + fe(2,3)*fe(2,3)
      be(3) = fe(3,1)*fe(3,1) + fe(3,2)*fe(3,2) + fe(3,3)*fe(3,3)
      be(4) = fe(1,1)*fe(2,1) + fe(1,2)*fe(2,2) + fe(1,3)*fe(2,3)
      be(5) = fe(1,1)*fe(3,1) + fe(1,2)*fe(3,2) + fe(1,3)*fe(3,3)
      be(6) = fe(2,1)*fe(3,1) + fe(2,2)*fe(3,2) + fe(2,3)*fe(3,3)

      I13e = (be(1) + be(2) + be(3))/3.d0

c...  calculate Cauchy stress
      do i = 1,ntens
        stress(i) = (detfe*k*(detfe-1.d0)*xi(i) + mu*J23e*(be(i) - I13e*xi(i)))/detfe
      enddo  

c...  tangent 
      vol = K*detfe*(2.d0*detfe - 1.d0) + 2.d0/3.d0*mu*J23e*I13e
      
      ddsdde(1,1)= (vol + 2.d0/3.d0*mu*J23e*be(1))/detf
      ddsdde(2,2)= (vol + 2.d0/3.d0*mu*J23e*be(2))/detf
      ddsdde(3,3)= (vol + 2.d0/3.d0*mu*J23e*be(3))/detf
      ddsdde(1,2)= (vol - 2.d0/3.d0*mu*J23e*(be(1)+be(2)))/detf
      ddsdde(1,3)= (vol - 2.d0/3.d0*mu*J23e*(be(1)+be(3)))/detf
      ddsdde(2,3)= (vol - 2.d0/3.d0*mu*J23e*(be(2)+be(3)))/detf
      ddsdde(1,4)= 1.d0/3.d0*mu*J23e*be(4)/detf
      ddsdde(2,4)= 1.d0/3.d0*mu*J23e*be(4)/detf
      ddsdde(3,4)=-2.d0/3.d0*mu*J23e*be(4)/detf
      ddsdde(1,5)= 1.d0/3.d0*mu*J23e*be(5)/detf
      ddsdde(2,5)=-2.d0/3.d0*mu*J23e*be(5)/detf
      ddsdde(3,5)= 1.d0/3.d0*mu*J23e*be(5)/detf
      ddsdde(1,6)=-2.d0/3.d0*mu*J23e*be(6)/detf
      ddsdde(2,6)= 1.d0/3.d0*mu*J23e*be(6)/detf
      ddsdde(3,6)= 1.d0/3.d0*mu*J23e*be(6)/detf
      ddsdde(4,4)= mu*J23e*(be(1)+be(2))/2.d0/detf
      ddsdde(5,5)= mu*J23e*(be(1)+be(3))/2.d0/detf
      ddsdde(6,6)= mu*J23e*(be(2)+be(3))/2.d0/detf
      ddsdde(4,5)= mu*J23e*be(6)/2.d0/detf
      ddsdde(4,6)= mu*J23e*be(5)/2.d0/detf
      ddsdde(5,6)= mu*J23e*be(4)/2.d0/detf

c...  use symmetry to fill in the rest
      do i=2, 6
          do j=1, i-1
              ddsdde(i,j)=ddsdde(j,i)
          end do
      end do      
     
c...  calculate strain energy 
      sse = (k*(detf - 1.d0)**2.d0 + 3.d0*mu*(I13e - 1.d0))/2.d0
     
      return
      end
