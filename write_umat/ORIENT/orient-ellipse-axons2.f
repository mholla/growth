c...  ------------------------------------------------------------------
      subroutine orient (T,noel,npt,layer,kspt,coords,basis,
     1 orname,nnodes,cnodes,jnnum)
c...  ------------------------------------------------------------------

      include 'aba_param.inc'

      character*80 orname 
      dimension T(3,3),coords(3),basis(3,3),cnodes(3,nnodes)
      dimension jnnum(nnodes)
      
      if      (orname(1:13).eq.'BRAIN-1-ORI-1') then
	  		call orient_ellipse(T,noel,npt,layer,kspt,coords,basis,orname,nnodes,cnodes,jnnum)
	  else if (orname(1:13).eq.'BRAIN-1-ORI-2') then
			call orient_axons(T,noel,npt,layer,kspt,coords,basis,orname,nnodes,cnodes,jnnum)
	  else 
			print *, "Please indicate a valid orientation name", orname
	  end if
	  
	  return 
	  end
	  
c...  ------------------------------------------------------------------
      subroutine orient_ellipse (T,noel,npt,layer,kspt,coords,basis,
     1 orname,nnodes,cnodes,jnnum)
c...  ------------------------------------------------------------------

      include 'aba_param.inc'

      character*80 orname
 
      dimension T(3,3),coords(3),basis(3,3),cnodes(3,nnodes)
      dimension jnnum(nnodes)
      
      real*8 AA,BB,a,b,abs_norm
      
      AA=2.d0
      BB=1.6d0
      b = sqrt(coords(1)**2.0d0*BB**2.0d0/AA**2.0d0+coords(2)**2.0d0)
      a = AA/BB*b
      
      if ((coords(1).eq.0.0d0).and.(coords(2).eq.0.0d0)) then
		  T(1,1) = 1.0d0
		  T(2,1) = 0.0d0
		  T(3,1) = 0.0d0
		  T(1,2) = 0.0d0
		  T(2,2) = 1.0d0
		  T(3,2) = 0.0d0
		  T(1,3) = 0.0d0
		  T(2,3) = 0.0d0
		  T(3,3) = 1.0d0 
      else
      	  abs_norm = sqrt((4.0d0*coords(1)**2.0d0/a**4.0d0)+(4.0d0*coords(2)**2.0d0/b**4.0d0))
		  T(1,1) = (2.0d0*coords(2)/b**2.0d0)/abs_norm
		  T(2,1) = (-2.0d0*coords(1)/a**2.0d0)/abs_norm
		  T(3,1) = 0.0d0
      
		  T(1,2) = (2.0d0*coords(1)/a**2.0d0)/abs_norm
		  T(2,2) = (2.0d0*coords(2)/b**2.0d0)/abs_norm
		  T(3,2) = 0.0d0
      
		  T(1,3) = 0.0d0
		  T(2,3) = 0.0d0
		  T(3,3) = 1.0d0
	  end if     
	  
      return
      end
      
c...  ------------------------------------------------------------------
      subroutine orient_axons (T,noel,npt,layer,kspt,coords,basis,
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