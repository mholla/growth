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
      subroutine orient (T,noel,npt,layer,kspt,coords,basis,orname,nnodes,cnodes,jnnum)
c...  ------------------------------------------------------------------
      include 'aba_param.inc'
C
      character*80 orname
C 
      dimension T(3,3),coords(3),basis(3,3),cnodes(3,nnodes),jnnum(nnodes)
C      
      T(1,1) = 1.0d0
      T(2,1) = 0.0d0
      T(3,1) = 0.0d0
      T(1,2) = 0.0d0
      T(2,2) = 1.0d0
      T(3,2) = 0.0d0
      T(1,3) = 0.0d0
      T(2,3) = 0.0d0
      T(3,3) = 1.0d0 
      
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

      call umat_area_morph_orient(stress,statev,ddsdde,sse,
     #                       time,dtime,coords,props,dfgrd1,
     #                       ntens,ndi,nshr,nstatv,nprops,
     #                       noel,npt,kstep,kinc)

      return
      end

c...  ------------------------------------------------------------------
      subroutine umat_area_morph_orient(stress,statev,ddsdde,sse,
     #                             time,dtime,coords,props,dfgrd1,
     #                             ntens,ndi,nshr,nstatv,nprops,
     #                             noel,npt,kstep,kinc)
c...  ------------------------------------------------------------------
c * " This UMAT is written for a Neohookean constitutive material.
c * " It implements morphogenetic area growth in the plane normal to xn0.
c * " It can ONLY be used with local material directions.
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
      real*8  lam, mu, alpha, xx(3)

c...  local variables      
      integer i,j,k
      real*8  norm, xn(3), xn0(3)
      real*8  theg, detf, detfe, lnJe
      real*8  fe(3,3), be(6), xi(6)
      real*8  b(6),finv(3,3),val(3),vec(3,3),V(3,3)
      
      data xi/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/

c...  initialize material parameters 
      lam   = props(1)              ! lame constant
      mu    = props(2)              ! shear modulus
      alpha = props(3)              ! growth rate
      xx(1) = props(4)              ! f0[1], unrotated
      xx(2) = props(5)              ! f0[2], unrotated
      xx(3) = props(6)              ! f0[3], unrotated

c...  calculate left stretch tensor v from F
      b(1) = dfgrd1(1,1)*dfgrd1(1,1) 
     #     + dfgrd1(1,2)*dfgrd1(1,2) 
     #     + dfgrd1(1,3)*dfgrd1(1,3)
      b(2) = dfgrd1(2,1)*dfgrd1(2,1) 
     #     + dfgrd1(2,2)*dfgrd1(2,2) 
     #     + dfgrd1(2,3)*dfgrd1(2,3)
      b(3) = dfgrd1(3,1)*dfgrd1(3,1) 
     #     + dfgrd1(3,2)*dfgrd1(3,2) 
     #     + dfgrd1(3,3)*dfgrd1(3,3)
      b(4) = dfgrd1(1,1)*dfgrd1(2,1) 
     #     + dfgrd1(1,2)*dfgrd1(2,2) 
     #     + dfgrd1(1,3)*dfgrd1(2,3)
      b(5) = dfgrd1(1,1)*dfgrd1(3,1) 
     #     + dfgrd1(1,2)*dfgrd1(3,2) 
     #     + dfgrd1(1,3)*dfgrd1(3,3)
      b(6) = dfgrd1(2,1)*dfgrd1(3,1) 
     #     + dfgrd1(2,2)*dfgrd1(3,2) 
     #     + dfgrd1(2,3)*dfgrd1(3,3)

      do i = 1,3
        do j = 1,3
          V(i,j) = 0.d0
        enddo
      enddo

      call sprind(b,val,vec,1,ndi,nshr)
      
      do k = 1,3
        do i = 1,3
          do j = 1,3
            V(i,j) = V(i,j) + dsqrt(val(k))*vec(k,i)*vec(k,j)
          enddo
        enddo
      enddo

c...  calculate deformed fiber in local rotated coordinate system
      norm = sqrt(xx(1)*xx(1) + xx(2)*xx(2) + xx(3)*xx(3))
      xn(1) = (V(1,1)*xx(1) + V(1,2)*xx(2) + V(1,3)*xx(3))/norm
      xn(2) = (V(2,1)*xx(1) + V(2,2)*xx(2) + V(2,3)*xx(3))/norm
      xn(3) = (V(3,1)*xx(1) + V(3,2)*xx(2) + V(3,3)*xx(3))/norm
     
c...  calculate F inverse
      detf = dfgrd1(1,1) * (dfgrd1(2,2)*dfgrd1(3,3)-dfgrd1(2,3)*dfgrd1(3,2))
     #     - dfgrd1(1,2) * (dfgrd1(2,1)*dfgrd1(3,3)-dfgrd1(2,3)*dfgrd1(3,1))
     #     + dfgrd1(1,3) * (dfgrd1(2,1)*dfgrd1(3,2)-dfgrd1(2,2)*dfgrd1(3,1))

      finv(1,1) = (dfgrd1(2,2)*dfgrd1(3,3) - dfgrd1(2,3)*dfgrd1(3,2))/detf
      finv(2,2) = (dfgrd1(1,1)*dfgrd1(3,3) - dfgrd1(1,3)*dfgrd1(3,1))/detf
      finv(3,3) = (dfgrd1(1,1)*dfgrd1(2,2) - dfgrd1(1,2)*dfgrd1(2,1))/detf
      finv(1,2) = (dfgrd1(1,3)*dfgrd1(3,2) - dfgrd1(1,2)*dfgrd1(3,3))/detf
      finv(1,3) = (dfgrd1(1,2)*dfgrd1(2,3) - dfgrd1(1,3)*dfgrd1(2,2))/detf
      finv(2,1) = (dfgrd1(2,3)*dfgrd1(3,1) - dfgrd1(2,1)*dfgrd1(3,3))/detf
      finv(2,3) = (dfgrd1(2,1)*dfgrd1(1,3) - dfgrd1(2,3)*dfgrd1(1,1))/detf
      finv(3,1) = (dfgrd1(3,2)*dfgrd1(2,1) - dfgrd1(3,1)*dfgrd1(2,2))/detf
      finv(3,2) = (dfgrd1(3,1)*dfgrd1(1,2) - dfgrd1(3,2)*dfgrd1(1,1))/detf

c...  calculate reference fiber in local rotated coordinate system
      xn0(1) = finv(1,1)*xn(1) + finv(1,2)*xn(2) + finv(1,3)*xn(3)
      xn0(2) = finv(2,1)*xn(1) + finv(2,2)*xn(2) + finv(2,3)*xn(3)
      xn0(3) = finv(3,1)*xn(1) + finv(3,2)*xn(2) + finv(3,3)*xn(3)

c...  update growth factor
      theg = 1.d0 + alpha*time(2)
      statev(1) = theg
      statev(2) = theg
      
c...  calculate elastic tensor Fe = 1\sqrt(theg)*F + (1-1/sqrt(theg))*n \otimes n_0
      fe(1,1) = 1.d0/sqrt(theg)*dfgrd1(1,1) + (1.d0 - 1.d0/sqrt(theg))*xn(1)*xn0(1)
      fe(1,2) = 1.d0/sqrt(theg)*dfgrd1(1,2) + (1.d0 - 1.d0/sqrt(theg))*xn(1)*xn0(2)
      fe(1,3) = 1.d0/sqrt(theg)*dfgrd1(1,3) + (1.d0 - 1.d0/sqrt(theg))*xn(1)*xn0(3)
      fe(2,1) = 1.d0/sqrt(theg)*dfgrd1(2,1) + (1.d0 - 1.d0/sqrt(theg))*xn(2)*xn0(1)
      fe(2,2) = 1.d0/sqrt(theg)*dfgrd1(2,2) + (1.d0 - 1.d0/sqrt(theg))*xn(2)*xn0(2)
      fe(2,3) = 1.d0/sqrt(theg)*dfgrd1(2,3) + (1.d0 - 1.d0/sqrt(theg))*xn(2)*xn0(3)
      fe(3,1) = 1.d0/sqrt(theg)*dfgrd1(3,1) + (1.d0 - 1.d0/sqrt(theg))*xn(3)*xn0(1)
      fe(3,2) = 1.d0/sqrt(theg)*dfgrd1(3,2) + (1.d0 - 1.d0/sqrt(theg))*xn(3)*xn0(2)
      fe(3,3) = 1.d0/sqrt(theg)*dfgrd1(3,3) + (1.d0 - 1.d0/sqrt(theg))*xn(3)*xn0(3)

c...  calculate determinant of deformation gradient
      detfe = +fe(1,1) * (fe(2,2)*fe(3,3)-fe(2,3)*fe(3,2))
     #        -fe(1,2) * (fe(2,1)*fe(3,3)-fe(2,3)*fe(3,1))
     #        +fe(1,3) * (fe(2,1)*fe(3,2)-fe(2,2)*fe(3,1))
     
      statev(3) = detfe
      lnJe = dlog(detfe)
     
c...  calculate elastic left cauchy-green deformation tensor be = fe * fe^t
      be(1) = fe(1,1)*fe(1,1) + fe(1,2)*fe(1,2) + fe(1,3)*fe(1,3)
      be(2) = fe(2,1)*fe(2,1) + fe(2,2)*fe(2,2) + fe(2,3)*fe(2,3)
      be(3) = fe(3,1)*fe(3,1) + fe(3,2)*fe(3,2) + fe(3,3)*fe(3,3)
      be(4) = fe(1,1)*fe(2,1) + fe(1,2)*fe(2,2) + fe(1,3)*fe(2,3)
      be(5) = fe(1,1)*fe(3,1) + fe(1,2)*fe(3,2) + fe(1,3)*fe(3,3)
      be(6) = fe(2,1)*fe(3,1) + fe(2,2)*fe(3,2) + fe(2,3)*fe(3,3)

c...  calculate Cauchy stress
      do i = 1,ntens
        stress(i) = ((lam*lnJe-mu)*xi(i) + mu*be(i))/detfe
      enddo   

c...  calculate elastic and geometric tangent
      ddsdde(1,1) = (lam + 2.d0*(lam*lnJe - mu))/detfe + 2.d0*stress(1)
      ddsdde(2,2) = (lam + 2.d0*(lam*lnJe - mu))/detfe + 2.d0*stress(2)
      ddsdde(3,3) = (lam + 2.d0*(lam*lnJe - mu))/detfe + 2.d0*stress(3)
      ddsdde(1,2) = (lam)/detfe
      ddsdde(1,3) = (lam)/detfe
      ddsdde(2,3) = (lam)/detfe
      ddsdde(1,4) = stress(4)
      ddsdde(2,4) = stress(4)
      ddsdde(3,4) = 0.d0
      ddsdde(4,4) = (lam*lnJe - mu)/detfe + (stress(1) + stress(2))/2.d0
      
      if (ntens.eq.6) then
        ddsdde(1,5) = stress(5)
        ddsdde(2,5) = 0.d0
        ddsdde(3,5) = stress(5)
        ddsdde(1,6) = 0.d0
        ddsdde(2,6) = stress(6)
        ddsdde(3,6) = stress(6)
        ddsdde(5,5) = (lam*lnJe - mu)/detfe + (stress(1) + stress(3))/2.d0
        ddsdde(6,6) = (lam*lnJe - mu)/detfe + (stress(2) + stress(3))/2.d0
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
      sse = (lam*lnJe**2.d0 + mu*(be(1)+be(2)+be(3) - 3.d0 - 2.d0*lnJe))/2.d0

      return
      end
