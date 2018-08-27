c...  ------------------------------------------------------------------
      subroutine umat_subcortex(stress,statev,ddsdde,sse,coords,time,dtime,props,dfgrd1)
c...  ------------------------------------------------------------------
      implicit none
	  
c...  variables passed in for information
	  real*8  props(6),dfgrd1(3,3),dtime,coords(3),time(2)
c...  variables to be defined
	  real*8 stress(6),ddsdde(6,6),statev(5),sse
c...  material properties
	  real*8 lam,mu,crMe,alpha,L,perturb,tf
c...  local variables
	  integer i,j,nitl
	  real*8 detf,fe(3,3),detfe
	  real*8 b(6),be(6),ce(3,3),ceinv(6),lnJe
	  real*8 theg,theg_n,sfac,dsfac,kfac,dkfac,res,dres
	  real*8 fac,cg(6)
	  real*8 center,delta
	  real*8 xi(6),xtol

	  
	  data xi/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/
	  xtol = 1.d-12
	  nitl = 0
	  
c...  ------------------------------------------------------------------
c...  initialize material parameters
c...  ------------------------------------------------------------------

      lam   = props(1)              ! lame constant
      mu    = props(2)              ! shear modulus
      alpha = props(3)				! 1/tau for negative growth
	  L     = props(4)				! length of rectangle (for perturbation)
	  perturb = props(5) 			! magnitude of perturbation
	  tf	= props(6)				! end point for perturbation
	  
	  center = L/2.
	  delta = L/40.
