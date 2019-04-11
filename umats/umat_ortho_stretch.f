c...  ------------------------------------------------------------------
      subroutine sdvini(statev,coords,nstatv,ncrds,noel,npt,layer,kspt)
c...  ------------------------------------------------------------------
      include 'aba_param.inc'

      dimension statev(nstatv)

      statev(1)=1.0d0    ! theg(1)
      statev(2)=1.0d0    ! theg(2)
      statev(3)=1.0d0    ! theg(3)
      statev(4)=1.0d0    ! the(1)
      statev(5)=1.0d0    ! the(2)
      statev(6)=1.0d0    ! the(3)
      statev(7)=1.0d0    ! Jg
      statev(8)=1.0d0    ! Je

      return
      end

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
     #ddsdde(ntens,ntens),ddsddt(ntens),drplde(ntens),
     #stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     #props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)

      call umat_ortho_stretch(stress,statev,ddsdde,sse,
     #                        time,dtime,coords,props,dfgrd1,
     #                        ntens,ndi,nshr,nstatv,nprops,
     #                        noel,npt,kstep,kinc)

      return
      end

c...  ------------------------------------------------------------------
      subroutine umat_ortho_stretch(stress,statev,ddsdde,sse,
     #                              time,dtime,coords,props,dfgrd1,
     #                              ntens,ndi,nshr,nstatv,nprops,
     #                              noel,npt,kstep,kinc)
c...  ------------------------------------------------------------------
c * " This UMAT is written for a Neohookean constitutive material.
c * " It implements orthotropic-stretch-driven growth.
c * " It is not suitable for use with local material directions.

c * " Written by Maria Holland in Aug 2012
c * " Updated by Maria Holland in Sep 2017
c...  ------------------------------------------------------------------

      implicit none

c...  variables to be defined
      real*8  stress(ntens),statev(nstatv),ddsdde(ntens,ntens),sse

c...  variables passed in for information
      real*8  time(2),dtime,coords(3),props(nprops),dfgrd1(3,3)
      integer ntens,ndi,nshr,nstatv,nprops,noel,npt,kstep,kinc

c...  material properties
      real*8  lam, tcr, mu, cr(3), max(3), alpha(3), cr_pos(3), cr_neg(3), phi_pos(3), phi_neg(3) 

c...  local variables      
      integer i,j,nitl
      real*8  xx0(3,3),xx(3,3),xxx(6,3),norm
      real*8  fginv(6),fe(3,3),detfe,be(6),lnJe
      real*8  detf,the(3),theg(3),theg_n(3)
      real*8  kg(3),dkg(3),phig(3),dphig(3),res(3),dres(3)
      real*8  fac1(3), fac2(3), cg_ij(6)
      real*8  xi(6),xtol

      data xi/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/
      xtol = 1.d-12

c...  initialize material parameters 
      lam       = props(1)             ! lame constant
      mu        = props(2)             ! shear modulus
      xx0(1,1)  = props(3)             ! r0[1]
      xx0(2,1)  = props(4)             ! r0[2]
      xx0(3,1)  = props(5)             ! r0[3]
      xx0(1,2)  = props(6)             ! t0[1]
      xx0(2,2)  = props(7)             ! t0[2]
      xx0(3,2)  = props(8)             ! t0[3]
      alpha(1)  = props(9)
      alpha(2)  = props(10)
      alpha(3)  = props(11)
      cr_pos(1) = props(12)
      cr_pos(2) = props(13)
      cr_pos(3) = props(14)
      cr_neg(1) = props(15)
      cr_neg(2) = props(16)
      cr_neg(3) = props(17)

c.... calculate undeformed vectors      
      norm = sqrt(xx0(1,1)*xx0(1,1) + xx0(2,1)*xx0(2,1) + xx0(3,1)*xx0(3,1))
      xx0(1,1) = xx0(1,1)/norm 
      xx0(2,1) = xx0(2,1)/norm 
      xx0(3,1) = xx0(3,1)/norm 
      
      norm = sqrt(xx0(1,2)*xx0(1,2) + xx0(2,2)*xx0(2,2) + xx0(3,2)*xx0(3,2))
      xx0(1,2) = xx0(1,2)/norm 
      xx0(2,2) = xx0(2,2)/norm 
      xx0(3,2) = xx0(3,2)/norm 
            
      xx0(1,3) = xx0(2,1)*xx0(3,2) - xx0(3,1)*xx0(2,2)
      xx0(2,3) = xx0(3,1)*xx0(1,2) - xx0(1,1)*xx0(3,2)
      xx0(3,3) = xx0(1,1)*xx0(2,2) - xx0(2,1)*xx0(1,2)

      detf = +dfgrd1(1,1)*(dfgrd1(2,2)*dfgrd1(3,3)-dfgrd1(2,3)*dfgrd1(3,2))
     #       -dfgrd1(1,2)*(dfgrd1(2,1)*dfgrd1(3,3)-dfgrd1(2,3)*dfgrd1(3,1))
     #       +dfgrd1(1,3)*(dfgrd1(2,1)*dfgrd1(3,2)-dfgrd1(2,2)*dfgrd1(3,1))

c...  growth in orthotropic directions
c...  ------------------------------------------------------------------
      do i = 1,3

c...    calculate deformed elastic normals
        xx(1,i) = dfgrd1(1,1)*xx0(1,i) + dfgrd1(1,2)*xx0(2,i) + dfgrd1(1,3)*xx0(3,i)
        xx(2,i) = dfgrd1(2,1)*xx0(1,i) + dfgrd1(2,2)*xx0(2,i) + dfgrd1(2,3)*xx0(3,i)
        xx(3,i) = dfgrd1(3,1)*xx0(1,i) + dfgrd1(3,2)*xx0(2,i) + dfgrd1(3,3)*xx0(3,i)

        xxx(1,i) = xx(1,i)*xx(1,i)
        xxx(2,i) = xx(2,i)*xx(2,i)
        xxx(3,i) = xx(3,i)*xx(3,i)
        xxx(4,i) = xx(1,i)*xx(2,i)
        xxx(5,i) = xx(1,i)*xx(3,i)
        xxx(6,i) = xx(2,i)*xx(3,i)

c       calculate stretch in i-th direction              
        the(i)  = sqrt(xx(1,i)*xx(1,i) + xx(2,i)*xx(2,i) + xx(3,i)*xx(3,i))
            
c       obtain state variable history
        theg_n(i) = statev(i)
        theg(i)   = theg_n(i)
              
c       local newton iteration 
c       ------------------------------------------------------------------
        nitl = 0

        phi_pos = the(i)/theg(i) - cr_pos(i) ! criterion for positive growth
        phi_neg = the(i)/theg(i) - cr_neg(i) ! criterion for negative growth
          
        if (phi_pos(i).gt.0.d0) then         ! positive growth
          phig(i) = phi_pos(i)
          tcr = cr_pos(i)
        else if (phi_neg(i).lt.0.d0) then    ! negative growth
          phig(i) = phi_neg(i)
          tcr = cr_neg(i)
        else                                 ! no growth
          phig(i) = 0.d0
        endif

        if (phig(i).eq.0.d0) then               ! no growth
          theg(i) = theg_n(i)
          fac1(i) = 0.d0  
        
        else                                 ! growth
200       continue
          nitl = nitl + 1

          kg(i)    = alpha(i)
          dkg(i)   = 0.d0
          phig(i)  = the(i)/theg(i) - cr(i)
          dphig(i) = -the(i)/theg(i)**2.d0
    
          res(i)   = theg(i) - theg_n(i) - kg(i)*phig(i)*dtime
          dres(i)  = 1.d0 - (kg(i)*dphig(i) + dkg(i)*phig(i))*dtime
                
          theg(i) = theg(i) - res(i) / dres(i)
            
          if ((nitl.lt.20).and.(dabs(res(i)).gt.xtol)) go to 200
          if (nitl.eq.20) print *, 'no local convergence in ',i, '! |r|=', res(i)
            
          fac1(i) = kg(i)*dtime/dres(i)/theg(i)/theg(i)/the(i)
          fac2(i) = 2.d0*mu/theg(i)/theg(i)
            
        endif  
c       ------------------------------------------------------------------
c       end local newton iteration 

c       update state variables
        statev(i) = theg(i)
        statev(i+3) = the(i)    
      
      enddo
c     ------------------------------------------------------------------
c...  end growth in orthotropic directions

c...  calculate elastic tensor Fe = F * Fg^{-1}
      fe(1,1) = xx(1,1)*xx0(1,1)/theg(1) + xx(1,2)*xx0(1,2)/theg(2) + xx(1,3)*xx0(1,3)/theg(3) 
      fe(1,2) = xx(1,1)*xx0(2,1)/theg(1) + xx(1,2)*xx0(2,2)/theg(2) + xx(1,3)*xx0(2,3)/theg(3) 
      fe(1,3) = xx(1,1)*xx0(3,1)/theg(1) + xx(1,2)*xx0(3,2)/theg(2) + xx(1,3)*xx0(3,3)/theg(3) 
      fe(2,1) = xx(2,1)*xx0(1,1)/theg(1) + xx(2,2)*xx0(1,2)/theg(2) + xx(2,3)*xx0(1,3)/theg(3) 
      fe(2,2) = xx(2,1)*xx0(2,1)/theg(1) + xx(2,2)*xx0(2,2)/theg(2) + xx(2,3)*xx0(2,3)/theg(3) 
      fe(2,3) = xx(2,1)*xx0(3,1)/theg(1) + xx(2,2)*xx0(3,2)/theg(2) + xx(2,3)*xx0(3,3)/theg(3) 
      fe(3,1) = xx(3,1)*xx0(1,1)/theg(1) + xx(3,2)*xx0(1,2)/theg(2) + xx(3,3)*xx0(1,3)/theg(3) 
      fe(3,2) = xx(3,1)*xx0(2,1)/theg(1) + xx(3,2)*xx0(2,2)/theg(2) + xx(3,3)*xx0(2,3)/theg(3) 
      fe(3,3) = xx(3,1)*xx0(3,1)/theg(1) + xx(3,2)*xx0(3,2)/theg(2) + xx(3,3)*xx0(3,3)/theg(3) 

c...  calculate determinant of elastic deformation gradient
      detfe = +fe(1,1)*(fe(2,2)*fe(3,3)-fe(2,3)*fe(3,2))
     #        -fe(1,2)*(fe(2,1)*fe(3,3)-fe(2,3)*fe(3,1))
     #        +fe(1,3)*(fe(2,1)*fe(3,2)-fe(2,2)*fe(3,1))
     
      statev(7) = theg(1)*theg(2)*theg(3)
      statev(8) = detfe
      lnJe = dlog(detfe)

c...  calculate left cauchy-green deformation tensor be = Fe * Fe^t
      be(1) = fe(1,1)*fe(1,1) + fe(1,2)*fe(1,2) + fe(1,3)*fe(1,3)
      be(2) = fe(2,1)*fe(2,1) + fe(2,2)*fe(2,2) + fe(2,3)*fe(2,3)
      be(3) = fe(3,1)*fe(3,1) + fe(3,2)*fe(3,2) + fe(3,3)*fe(3,3)
      be(4) = fe(1,1)*fe(2,1) + fe(1,2)*fe(2,2) + fe(1,3)*fe(2,3)
      be(5) = fe(1,1)*fe(3,1) + fe(1,2)*fe(3,2) + fe(1,3)*fe(3,3)
      be(6) = fe(2,1)*fe(3,1) + fe(2,2)*fe(3,2) + fe(2,3)*fe(3,3)

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
      
c...  compile tangent
      cg_ij(1) = (lam*lnJe - lam - mu)*xi(1) + mu*be(1)
      cg_ij(2) = (lam*lnJe - lam - mu)*xi(2) + mu*be(2)
      cg_ij(3) = (lam*lnJe - lam - mu)*xi(3) + mu*be(3)
      cg_ij(4) = (lam*lnJe - lam - mu)*xi(4) + mu*be(4)
      cg_ij(5) = (lam*lnJe - lam - mu)*xi(5) + mu*be(5)
      cg_ij(6) = (lam*lnJe - lam - mu)*xi(6) + mu*be(6)

      do i = 1,ntens
        do j = 1,ntens
          ddsdde(i,j) = ddsdde(i,j) 
     &                + fac1(1)*( cg_ij(i) - fac2(i)*xxx(i,1) )*xxx(j,1)/detfe
     &                + fac1(2)*( cg_ij(i) - fac2(i)*xxx(i,2) )*xxx(j,2)/detfe
     &                + fac1(3)*( cg_ij(i) - fac2(i)*xxx(i,3) )*xxx(j,3)/detfe
        enddo
      enddo

c...  calculate strain energy
      sse = (lam*lnJe**2.d0 + mu*(be(1)+be(2)+be(3) - 3.d0 - 2.d0*lnJe))/2.d0

      return
      end
