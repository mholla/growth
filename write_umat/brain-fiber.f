c...  ------------------------------------------------------------------
      subroutine umat_subcortex(stress,statev,ddsdde,sse,coords,time,dtime,props,dfgrd1,ntens,ndi,nstatev)
c...  ------------------------------------------------------------------
      implicit none
      
c...  variables passed in for information
      real*8  props(10),dfgrd1(3,3),dtime,coords(3),time(2)
      integer ntens,ndi,nstatev
c...  variables to be defined
      real*8 stress(6),ddsdde(6,6),statev(5),sse
c...  material properties
      real*8 lam,mu,alpha,tcr,L,perturb,tf
c...  local variables
      integer i,j,nitl,noel
      real*8 xn(3),nn(6),xnt(3),nnt(6),norm,ftn(3)
      real*8 detf,fginv(6),fe(3,3),detfe,lnJe,be(6)
      real*8 the,theg,theg_n,phi,fac
      real*8 cg_ij(6),xi(6),xtol,center,delta
      
      data xi/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/
      xtol = 1.d-12
      
c...  ------------------------------------------------------------------
c...  initialize material parameters
c...  ------------------------------------------------------------------
      lam   = props(1)              ! lame constant
      mu    = props(2)              ! shear modulus
      alpha = props(3)
      tcr   = props(4)
      xn(1) = props(5)              ! n0[1]
      xn(2) = props(6)              ! n0[2]
      xn(3) = props(7)              ! n0[3]
            
      L       = props(8)        
      perturb = props(9)            
      tf      = props(10)
      center  = L/2.d0
      delta   = L/4.d1

