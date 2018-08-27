c...  ------------------------------------------------------------------
      subroutine orient (T,noel,npt,layer,kspt,coords,basis,
     1 orname,nnodes,cnodes,jnnum)
c...  ------------------------------------------------------------------

      include 'aba_param.inc'

      character*80 orname
 
      dimension T(3,3),coords(3),basis(3,3),cnodes(3,nnodes)
      dimension jnnum(nnodes)
      
      real*8 zero, one, two, pi
      real*8 HH,LL,AA,BB,r,a,b
      real*8 theta,rr(3)
      
      zero = 0.d0
      one = 1.d0
      two = 2.d0
      pi = 4.*atan(1.)
      
      HH = 5.d-1
      LL = 12.d-1
      
      AA = LL
      BB = HH*two

      r = sqrt((coords(1)/AA)**two + (coords(2)/BB)**two)    
	  b  = r**two * BB*two
	  a  = AA/BB*b

c...  This subroutine sets axons tangentially at the cortex, radially at the core, with a smooth transition
	  if (coords(1).le.zero) then
		  theta = - pi/two/HH*r
      else
		  theta =   pi/two/HH*r
	  end if 

c...  radial line	  
	  rr(1) = two*coords(1)/a**two
      rr(2) = two*coords(2)/b**two
	  norm = sqrt(rr(1)**two + rr(2)**two)      
            
	  T(1,1) = rr(1)*cos(theta) - rr(2)*sin(theta)
	  T(2,1) = rr(1)*sin(theta) + rr(2)*cos(theta)
	  T(3,1) = zero
  
	  T(1,2) = zero
	  T(2,2) = zero
	  T(3,2) = one

      return
      end