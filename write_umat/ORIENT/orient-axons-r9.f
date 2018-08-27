c...  ------------------------------------------------------------------
      subroutine orient (T,noel,npt,layer,kspt,coords,basis,
     1 orname,nnodes,cnodes,jnnum)
c...  ------------------------------------------------------------------

      include 'aba_param.inc'

      character*80 orname
 
      dimension T(3,3),coords(3),basis(3,3),cnodes(3,nnodes)
      dimension jnnum(nnodes)
      
      real*8 zero, one, two, pi
      real*8 HH,LL
      real*8 XX(3),theta
      
      zero = 0.d0
      one = 1.d0
      two = 2.d0
      pi = 4.*atan(1.)
      
      HH = 1.d0
      LL = 2.d0
      
c... This subroutine sets axons tangentially at the cortex, radially at the core, with a smooth transition
      
c...  horizontal line      
      XX(1) = one
      XX(2) = zero
      XX(3) = zero
      
      if (coords(1).le.LL/two) then
		  theta = - pi/two/HH*coords(2)
      else
		  theta =   pi/two/HH*coords(2)
	  end if 
	  
	  T(1,1) = XX(1)*cos(theta) - XX(2)*sin(theta)
	  T(2,1) = XX(1)*sin(theta) + XX(2)*cos(theta)
	  T(3,1) = XX(3)
  
	  T(1,2) = zero
	  T(2,2) = zero
	  T(3,2) = one

      return
      end