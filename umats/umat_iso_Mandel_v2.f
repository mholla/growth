
c...  ------------------------------------------------------------------
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     #rpl,ddsddt,drplde,drpldt,
     #stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     #ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     #celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
c...  ------------------------------------------------------------------
      include 'aba_param.inc'

      character*80 cmname
      dimension stress(ntens),statev(nstatv),
     # ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     # stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     # props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      call umat_iso_Mandel(stress,statev,ddsdde,sse,
     #                     time,dtime,coords,props,dfgrd1,
     #                     ntens,ndi,nshr,nstatv,nprops,
     #                     noel,npt,kstep,kinc)


      return
      end

c...  ------------------------------------------------------------------
      subroutine umat_iso_Mandel(stress,statev,ddsdde,sse,
     #                           time,dtime,coords,props,dfgrd1,
     #                           ntens,ndi,nshr,nstatv,nprops,
     #                           noel,npt,kstep,kinc)

c...  ------------------------------------------------------------------
c * " This UMAT is written for a Neohookean constitutive material.
c * " It implements Mandel-stress-driven isotropic growth.

c * " Written by Maria Holland in Aug 2012
c * " Updated by Maria Holland in Sep 2017
c * Updated by Karan Taneja in February 2024
c...  ------------------------------------------------------------------

      implicit none
      
c...  variables to be defined
      real*8  stress(ntens),statev(nstatv),ddsdde(ntens,ntens),sse

c...  variables passed in for information
      real*8  time(2),dtime,coords(3),props(nprops),dfgrd1(3,3)
      integer ntens,ndi,nshr,nstatv,nprops,noel,npt,kstep,kinc

c...  material properties
      real*8  lam, mu, alpha, cr_pos, cr_neg

c...  local variables
      integer i,j,nitl
      real*8  the, theg, theg_n, detf
      real*8  fe(3,3), be(6), detfe, lnJe
      real*8  phig, dphig, phi_pos, phi_neg, kg, dkg, res, dres
      real*8  trMe, trMe_cr, fac, cg_ij(6), cg_kl(6)
      real*8  xi(6), xtol


      ! Parameters
      !
      real*8 zero,one,two,half,three,third,nine
      parameter(zero=0.d0,one=1.d0,two=2.d0,half=0.5d0,three=3.d0,
     +     third=1.d0/3.d0,nine=9.d0)

      data xi/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/
      xtol = 1.d-12

c...  initialize material parameters 
      lam     = props(1)              ! lame constant
      mu      = props(2)              ! shear modulus
      alpha   = props(3)              ! growth rate
      cr_pos  = props(4)              ! critical stress for positive growth
      cr_neg  = props(5)              ! critical stress for negative growth
      
c...  calculate determinant of deformation gradient
      detf = +dfgrd1(1,1)*(dfgrd1(2,2)*dfgrd1(3,3)-dfgrd1(2,3)*dfgrd1(3,2))
     #       -dfgrd1(1,2)*(dfgrd1(2,1)*dfgrd1(3,3)-dfgrd1(2,3)*dfgrd1(3,1))
     #       +dfgrd1(1,3)*(dfgrd1(2,1)*dfgrd1(3,2)-dfgrd1(2,2)*dfgrd1(3,1))

      the = detf

c...  obtain state variable history
! KT: Changed this update the growth variable 
      if((kinc.le.1).and.(kstep.eq.1)) then

         
            theg_n = one
   
      else

            theg_n = statev(1)
         !
      endif
!     theg_n = statev(1)
      theg   = theg_n

c...  local newton iteration 
c...  ------------------------------------------------------------------
      nitl = 0
      
 200    continue      
     
c...    elastic tensor Fe = F * Fg^{-1}
        fe(1,1)=dfgrd1(1,1)/theg
        fe(1,2)=dfgrd1(1,2)/theg
        fe(1,3)=dfgrd1(1,3)/theg
        fe(2,1)=dfgrd1(2,1)/theg
        fe(2,2)=dfgrd1(2,2)/theg
        fe(2,3)=dfgrd1(2,3)/theg
        fe(3,1)=dfgrd1(3,1)/theg
        fe(3,2)=dfgrd1(3,2)/theg
        fe(3,3)=dfgrd1(3,3)/theg
 
        detfe = +fe(1,1)*(fe(2,2)*fe(3,3)-fe(2,3)*fe(3,2))
     #           -fe(1,2)*(fe(2,1)*fe(3,3)-fe(2,3)*fe(3,1))
     #           +fe(1,3)*(fe(2,1)*fe(3,2)-fe(2,2)*fe(3,1))
        lnJe = log(detfe)
        
c...    left cauchy-green deformation tensor be = Fe * Fe^t
        be(1) = fe(1,1)*fe(1,1) + fe(1,2)*fe(1,2) + fe(1,3)*fe(1,3)
        be(2) = fe(2,1)*fe(2,1) + fe(2,2)*fe(2,2) + fe(2,3)*fe(2,3)
        be(3) = fe(3,1)*fe(3,1) + fe(3,2)*fe(3,2) + fe(3,3)*fe(3,3)
        be(4) = fe(1,1)*fe(2,1) + fe(1,2)*fe(2,2) + fe(1,3)*fe(2,3)
        be(5) = fe(1,1)*fe(3,1) + fe(1,2)*fe(3,2) + fe(1,3)*fe(3,3)
        be(6) = fe(2,1)*fe(3,1) + fe(2,2)*fe(3,2) + fe(2,3)*fe(3,3)
             
c...    growth criterion
        trMe  = 3*(lam*lnJe - mu) + mu*(be(1)+be(2)+be(3))
        phi_pos =  trMe - cr_pos
        phi_neg =  trMe - cr_neg
        
c...    growth 
        if (phi_pos.gt.0.d0) then      ! positive growth
          trMe_cr = cr_pos
        elseif (phi_neg.lt.0.d0) then  ! negative growth
          trMe_cr = cr_neg
        else                           ! no growth
          trMe_cr = 1.234d0
        endif

        if (trMe_cr.eq.1.234d0) then   ! no growth
          res   = 0.d0  
          dres  = 1.d0
          fac   = 0.d0
        else                           ! growth
          nitl = nitl + 1
          phig = trMe - trMe_cr
          dphig =  -1.d0/3.d0/theg*(2.d0*mu*(be(1)+be(2)+be(3)) + 9.d0*lam)
          kg  = alpha
          dkg = 0.d0
          
          res   = theg - theg_n - kg*phig*dtime
          dres  = 1.d0 - (dkg*phig + kg*dphig)*dtime
          theg = theg - res / dres
          
          fac   = kg*dtime/detfe/dres      
  
c...      check for convergence      
          if ((nitl.lt.50).and.(dabs(res).gt.xtol)) go to 200
          if (nitl.eq.50) print *, '>no convergence!!!! |r|=', dabs(res)

        endif  
c...  ------------------------------------------------------------------
c...  end local newton iteration 

c...  update state variables 
      statev(1) = theg
      statev(2) = trMe
      statev(3) = phi_pos
      statev(4) = phi_neg
      statev(5) = trMe_cr

            
c...  calculate Cauchy stress
      do i = 1,ntens
        stress(i) = ((lam*lnJe-mu)*xi(i) + mu*be(i))/detfe
      enddo      
      
c...  calculate elastic and geometric tangent
      ddsdde(1,1) = (lam - 2.d0*(lam*lnJe - mu))/detfe + 2.d0*stress(1)
      ddsdde(2,2) = (lam - 2.d0*(lam*lnJe - mu))/detfe + 2.d0*stress(2)
      ddsdde(3,3) = (lam - 2.d0*(lam*lnJe - mu))/detfe + 2.d0*stress(3)
      ddsdde(1,2) = (lam)/detfe
      ddsdde(1,3) = (lam)/detfe
      ddsdde(2,3) = (lam)/detfe
      ddsdde(1,4) = stress(4)
      ddsdde(2,4) = stress(4)
      ddsdde(3,4) = 0.d0
      ddsdde(4,4) = -(lam*lnJe - mu)/detfe + (stress(1) + stress(2))/2.d0
      
      if (ntens.eq.6) then
        ddsdde(1,5) = stress(5)
        ddsdde(2,5) = 0.d0
        ddsdde(3,5) = stress(5)
        ddsdde(1,6) = 0.d0
        ddsdde(2,6) = stress(6)
        ddsdde(3,6) = stress(6)
        ddsdde(5,5) = -(lam*lnJe - mu)/detfe + (stress(1) + stress(3))/2.d0
        ddsdde(6,6) = -(lam*lnJe - mu)/detfe + (stress(2) + stress(3))/2.d0
        ddsdde(4,5) = stress(6)/2.d0
        ddsdde(4,6) = stress(5)/2.d0
        ddsdde(5,6) = stress(4)/2.d0
      endif
    
c...  use symmetry to fill in the rest
      do i = 2,ntens
        do j = 1,i-1
          ddsdde(i,j) = ddsdde(j,i)
        enddo
      enddo    

c...  calculate growth tangent
      cg_ij(1) = 3.d0*(lam*lnJe - lam - mu)*xi(1) + mu*be(1)
      cg_ij(2) = 3.d0*(lam*lnJe - lam - mu)*xi(2) + mu*be(2)
      cg_ij(3) = 3.d0*(lam*lnJe - lam - mu)*xi(3) + mu*be(3)
      cg_ij(4) = 3.d0*(lam*lnJe - lam - mu)*xi(4) + mu*be(4)
      cg_ij(5) = 3.d0*(lam*lnJe - lam - mu)*xi(5) + mu*be(5)
      cg_ij(6) = 3.d0*(lam*lnJe - lam - mu)*xi(6) + mu*be(6)

      cg_kl(1) = 3.d0*lam*xi(1) + 2.d0*mu*be(1)
      cg_kl(2) = 3.d0*lam*xi(2) + 2.d0*mu*be(2)
      cg_kl(3) = 3.d0*lam*xi(3) + 2.d0*mu*be(3)
      cg_kl(4) = 3.d0*lam*xi(4) + 2.d0*mu*be(4)
      cg_kl(5) = 3.d0*lam*xi(5) + 2.d0*mu*be(5)
      cg_kl(6) = 3.d0*lam*xi(6) + 2.d0*mu*be(6)

c...  compile tangent
      do i = 1,ntens
        do j = 1,ntens
          ddsdde(i,j) = ddsdde(i,j) + fac*cg_ij(i)*cg_kl(j)
        enddo
      enddo

c...  calculate strain energy
      sse = (lam*lnJe**2.d0 + mu*(be(1)+be(2)+be(3) - 3.d0 - 2.d0*lnJe))/2.d0

      statev(6) = sse

      return
      end
