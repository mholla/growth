      subroutine orient (T,noel,npt,layer,kspt,coords,basis,
     1 orname,nnodes,cnodes,jnnum)
C
      include 'aba_param.inc'
C
      character*80 orname
C 
      dimension T(3,3),coords(3),basis(3,3),cnodes(3,nnodes)
      dimension jnnum(nnodes)
	  
	  T(1,1) = cnodes(1,1)-cnodes(1,4)
	  T(2,1) = cnodes(2,1)-cnodes(2,4)
	  T(3,1) = cnodes(3,1)-cnodes(3,4)
	  
	  T(1,2) = cnodes(1,8)-cnodes(1,4)
	  T(2,2) = cnodes(2,8)-cnodes(2,4)
	  T(3,2) = cnodes(3,8)-cnodes(3,4)
C
      return
      end