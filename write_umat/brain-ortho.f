c...  ------------------------------------------------------------------
      subroutine umat_subcortex(stress,statev,ddsdde,sse,coords,time,dtime,props,dfgrd1,noel)
c...  ------------------------------------------------------------------
      implicit none

c...  variables passed in for information
      real*8  props(15),dfgrd1(3,3),dtime,coords(3),time(2)
c...  variables to be defined
      real*8  stress(6),ddsdde(6,6),statev(7),sse
c...  material properties
      real*8  lam,mu,cr(3),max(3),gam,alpha,L,perturb,tf
c...  local variables      
      integer i,j,nitl,noel
      real*8  xx(3,3),xxt(3,3),xox(6,3),xoxt(6,3),norm
      real*8  detf,finv(3,3),b(6)
      real*8  fginv(6),fe(3,3),detfe,be(6),lnJe
      real*8  the(3),theg(3),theg_n(3)
      real*8  phi(3),k(3),s(3),dk(3),ds(3),R(3),dR(3),fac(3)
      real*8  cg_ij(6,3),xi(6),xtol,center,delta
      
      data xi/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/
      xtol = 1.d-12

c...  initialize material parameters from ud-field
      lam      = props(1)           ! lame constant
      mu      = props(2)            ! shear modulus
      cr(1)   = props(3)
      cr(2)   = props(4)
      cr(3)   = props(5)
      alpha   = props(6)            ! 1/tau
      xx(1,1) = props(7)            ! r0[1]
      xx(2,1) = props(8)            ! r0[2]
      xx(3,1) = props(9)            ! r0[3]
      xx(1,2) = props(10)           ! t0[1]
      xx(2,2) = props(11)           ! t0[2]
      xx(3,2) = props(12)           ! t0[3]
            
      L       = props(13)        
      perturb = props(14)            
      tf      = props(15)
      center  = L/2.d0
      delta   = L/4.d1

