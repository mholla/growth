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
      subroutine umat_neohooke_Abaqus(stress,statev,ddsdde,sse,
     #                                time,dtime,coords,props,dfgrd1,
     #                                ntens,ndi,nshr,nstatv,nprops,
     #                                noel,npt,kstep,kinc)
c...  ------------------------------------------------------------------
c * " This UMAT is written for the Abaqus Neohookean constitutive material.
c * " It should yield the same results as the built-in model

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
      real*8  lam,mu, k

c...  local variables      
      integer i,j
      real*8  detf, I13, J53, vol
      real*8  b(6), xi(6)
      
      data xi/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/

c...  initialize material parameters
      lam   = props(1)              ! lame constant
      mu    = props(2)              ! shear modulus
      k = (3.d0*lam + 2.d0*mu)/3.d0
      
c...  calculate determinant of deformation gradient
      detf=+dfgrd1(1,1) * (dfgrd1(2,2)*dfgrd1(3,3)-dfgrd1(2,3)*dfgrd1(3,2))
     #     -dfgrd1(1,2) * (dfgrd1(2,1)*dfgrd1(3,3)-dfgrd1(2,3)*dfgrd1(3,1))
     #     +dfgrd1(1,3) * (dfgrd1(2,1)*dfgrd1(3,2)-dfgrd1(2,2)*dfgrd1(3,1))

      J53 = detf**(-5.d0/3.d0)
     
c...  left cauchy-green deformation tensor b = f * f^t
      b(1) = dfgrd1(1,1)*dfgrd1(1,1) + dfgrd1(1,2)*dfgrd1(1,2) + dfgrd1(1,3)*dfgrd1(1,3)
      b(2) = dfgrd1(2,1)*dfgrd1(2,1) + dfgrd1(2,2)*dfgrd1(2,2) + dfgrd1(2,3)*dfgrd1(2,3)
      b(3) = dfgrd1(3,1)*dfgrd1(3,1) + dfgrd1(3,2)*dfgrd1(3,2) + dfgrd1(3,3)*dfgrd1(3,3)
      b(4) = dfgrd1(1,1)*dfgrd1(2,1) + dfgrd1(1,2)*dfgrd1(2,2) + dfgrd1(1,3)*dfgrd1(2,3)
      b(5) = dfgrd1(1,1)*dfgrd1(3,1) + dfgrd1(1,2)*dfgrd1(3,2) + dfgrd1(1,3)*dfgrd1(3,3)
      b(6) = dfgrd1(2,1)*dfgrd1(3,1) + dfgrd1(2,2)*dfgrd1(3,2) + dfgrd1(2,3)*dfgrd1(3,3)

      I13 = (b(1) + b(2) + b(3))/3.d0

c...  Cauchy stress
      do i = 1,ntens
        stress(i) = k*(detf-1.d0)*xi(i) + mu*J53*(b(i) - I13*xi(i))
      enddo

c...  tangent 
      vol = K*(2.d0*detf - 1.d0) + 2.d0/3.d0*mu*J53*I13
      
      ddsdde(1,1)= vol + 2.d0/3.d0*mu*J53*b(1)
      ddsdde(2,2)= vol + 2.d0/3.d0*mu*J53*b(2)
      ddsdde(3,3)= vol + 2.d0/3.d0*mu*J53*b(3)
      ddsdde(1,2)= vol - 2.d0/3.d0*mu*J53*(b(1)+b(2))
      ddsdde(1,3)= vol - 2.d0/3.d0*mu*J53*(b(1)+b(3))
      ddsdde(2,3)= vol - 2.d0/3.d0*mu*J53*(b(2)+b(3))
      ddsdde(1,4)= 1.d0/3.d0*mu*J53*b(4)
      ddsdde(2,4)= 1.d0/3.d0*mu*J53*b(4)
      ddsdde(3,4)=-2.d0/3.d0*mu*J53*b(4)
      ddsdde(4,4)= mu*J53*(b(1)+b(2))/2.d0

      if (ntens.eq.6) then
        ddsdde(1,5)= 1.d0/3.d0*mu*J53*b(5)
        ddsdde(2,5)=-2.d0/3.d0*mu*J53*b(5)
        ddsdde(3,5)= 1.d0/3.d0*mu*J53*b(5)
        ddsdde(1,6)=-2.d0/3.d0*mu*J53*b(6)
        ddsdde(2,6)= 1.d0/3.d0*mu*J53*b(6)
        ddsdde(3,6)= 1.d0/3.d0*mu*J53*b(6)
        ddsdde(5,5)= mu*J53*(b(1)+b(3))/2.d0
        ddsdde(6,6)= mu*J53*(b(2)+b(3))/2.d0
        ddsdde(4,5)= mu*J53*b(6)/2.d0
        ddsdde(4,6)= mu*J53*b(5)/2.d0
        ddsdde(5,6)= mu*J53*b(4)/2.d0
      endif

c...  Use symmetry to fill in the rest
      do i=2, 6
        do j=1, i-1
          ddsdde(i,j)=ddsdde(j,i)
        end do
      end do      

c...  calculate strain energy 
      sse = (k*(detf - 1.d0)**2.d0 + 3.d0*mu*(I13 - 1.d0))/2.d0
     
      return
      end
