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
      
      HH = 5.d-1
      LL = 12.d-1
      
c...  This subroutine sets axons in a concave down orientation
      
	  if (coords(1).le.zero) then
		  cc(1) = -LL/two
		  cc(2) = zero
		  cc(3) = zero     
      else
		  cc(1) = LL/two
		  cc(2) = zero
		  cc(3) = zero
	  end if 
	  
	  rr(1) = coords(1) - cc(1)
	  rr(2) = coords(2) - cc(2)
	  norm = sqrt(rr(1)**two + rr(2)**two)

      if ((rr(1).eq.zero).and.(rr(2).eq.zero)) then
		  T(1,1) = one
		  T(2,1) = zero
		  T(3,1) = zero
		  T(1,2) = zero
		  T(2,2) = one
		  T(3,2) = zero
		  T(1,3) = zero
		  T(2,3) = zero
		  T(3,3) = one
      else      	  
		  T(1,1) =  rr(2)/norm
		  T(2,1) = -rr(1)/norm
		  T(3,1) =  zero
      
		  T(1,2) = rr(1)/norm
		  T(2,2) = rr(2)/norm
		  T(3,2) = zero
      
		  T(1,3) = zero
		  T(2,3) = zero
		  T(3,3) = one
	  end if     

      return
      end