c...  ------------------------------------------------------------------
      subroutine orient (T,noel,npt,layer,kspt,coords,basis,
     1 orname,nnodes,cnodes,jnnum)
c...  ------------------------------------------------------------------

      include 'aba_param.inc'

      character*80 orname
 
      dimension T(3,3),coords(3),basis(3,3),cnodes(3,nnodes)
      dimension jnnum(nnodes)
      
      real*8 zero, one, two
      real*8 HH,LL,norm
      real*8 cc(3),rr(3)
      
      zero = 0.d0
      one = 1.d0
      two = 2.d0
      
      HH = 1.d0
      LL = 2.d0
      
c... This subroutine sets all axons tangentially
      
	  T(1,1) = one
	  T(2,1) = zero
	  T(3,1) = zero
	  T(1,2) = zero
	  T(2,2) = one
	  T(3,2) = zero
	  T(1,3) = zero
	  T(2,3) = zero
	  T(3,3) = one
      
      return
      end