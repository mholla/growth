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
      
      AA=1.2d0
      BB=1.0d0
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
      
      real*8 zero, one, two, three, five
      real*8 HH,LL,AA,BB,r,a,b,norm
      real*8 rr(3)
      
      zero = 0.d0
      one = 1.d0
      two = 2.d0
      three = 3.d0
      five = 5.d0
      
      HH = 5.d-1
      LL = 12.d-1
      
      AA = LL
      BB = HH*two
      
      r = sqrt((coords(1)/AA)**two + (coords(2)/BB)**two)      
	  b  = r**two * BB*two
	  a  = AA/BB*b
	  
	  rr(1) = two*coords(1)/a**two
      rr(2) = two*coords(2)/b**two
	  norm = sqrt(rr(1)**two + rr(2)**two)
      
c... This subroutine sets axons tangentially near the cortex, and radially at the core
      
	  if (r.ge.two/five) then 	! tangential
		  T(1,1) =  rr(2)/norm
		  T(2,1) = -rr(1)/norm
		  T(3,1) =  zero
      
		  T(1,2) =  rr(1)/norm
		  T(2,2) =  rr(2)/norm
		  T(3,2) =  zero
      
      else 						 ! radial
		  T(1,1) =  rr(1)/norm
		  T(2,1) =  rr(2)/norm
		  T(3,1) =  zero
		  
		  T(1,2) =  rr(2)/norm
		  T(2,2) = -rr(1)/norm
		  T(3,2) =  zero

	  end if 
	  	  T(1,3) =  zero
		  T(2,3) =  zero
		  T(3,3) =  one

	  
      return
      end