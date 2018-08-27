c...  ------------------------------------------------------------------
      subroutine umat_subcortex(stress,statev,ddsdde,sse,dtime,props,dfgrd1,ntens,ndi,nstatv,nprops,noel,npt,kstep,kinc,coords)
c...  ------------------------------------------------------------------
      implicit none
      
c...  variables passed in for information
      integer ntens,ndi,nstatv,nprops,noel,npt,kstep,kinc
      real*8  props(nprops),dfgrd1(3,3),dtime,time(2),coords(3)
c...  variables to be defined
      real*8  stress(ntens),ddsdde(ntens,ntens),statev(nstatv),sse
c...  material properties
      real*8 lam,mu,alpha,cr_pos,cr_neg,L,perturb,tf
c...  local variables
      integer i,j,nitl
      real*8 xn(3),nn(6),xnt(3),nnt(6),norm,ftn(3)
      real*8 detf,fginv(6),fe(3,3),detfe,lnJe,be(6)
      real*8 the,theg,theg_n
      real*8 phi,phi_pos,phi_neg,crit
      real*8 kfac,sfac,dkfac,dsfac,res,dres,fac
      real*8 cg_ij(6),xi(6),xtol
      real*8 center,delta
      real*8 zero, one, two, HH, LL, cc(3), rr(3)

      data xi/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/
      xtol = 1.d-12
      zero = 0.d0
      one = 1.d0
      two = 2.d0
      HH = 1.d0
      LL = 2.d0      
      
c...  ------------------------------------------------------------------
c...  initialize material parameters
c...  ------------------------------------------------------------------
      lam    = props(1)             ! lame constant
      mu     = props(2)             ! shear modulus
      alpha  = props(3)
      cr_pos = props(4)
      cr_neg = props(5)
      
      L       = props(6)        
      perturb = props(7)            
      tf      = props(8)
      center  = L/2.d0
      delta   = L/4.d1

