c...  ------------------------------------------------------------------
      subroutine sdvini(statev,coords,nstatv,ncrds,noel,npt,layer,kspt)
c...  ------------------------------------------------------------------
      include 'aba_param.inc'

      dimension statev(nstatv)

      statev(1)=1.0d0    ! the_g
      statev(2)=1.0d0    ! the_e
      statev(3)=1.0d0    ! the
      statev(4)=1.0d0    ! xn(1)
      statev(5)=1.0d0    ! xn(2)
      statev(6)=1.0d0    ! xn(3)

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

      call umat_area_stretch(stress,statev,ddsdde,sse,
     #                       time,dtime,coords,props,dfgrd1,
     #                       ntens,ndi,nshr,nstatv,nprops,
     #                       noel,npt,kstep,kinc)

      return
      end

c...  ------------------------------------------------------------------
      subroutine umat_area_stretch(stress,statev,ddsdde,sse,
     #                             time,dtime,coords,props,dfgrd1,
     #                             ntens,ndi,nshr,nstatv,nprops,
     #                             noel,npt,kstep,kinc)
c...  ------------------------------------------------------------------
c * " This UMAT is written for a Neohookean constitutive material.
c * " It implements area-stretch-driven area growth in the plane normal to xn0.
c * " It is not suitable for use with local material directions.

c * " Written by Maria Holland in Aug 2012
c * " Updated by Maria Holland in Sep 2017
c...  ------------------------------------------------------------------

      implicit none

c...  variables to be defined
      real*8  stress(ntens), ddsdde(ntens,ntens), statev(nstatv), sse

c...  variables passed in for information
      real*8  time(2), dtime, coords(3), props(nprops), dfgrd1(3,3)
      integer ntens, ndi, nshr, nstatv, nprops, noel, npt, kstep, kinc

c...  material properties
      real*8  lam, mu, xn0(3), alpha, tcr, tmax, gam

c...  local variables
      integer i, j, nitl
      real*8  norm, xn(3), finv(3,3), ftn0(3)
      real*8  the, theg, theg_n, detf
      real*8  fe(3,3), be(6), detfe, lnJe
      real*8  kg, phig, dkg, dphig, res, dres
      real*8  fac, cg_ij(6), cg_kl(6)
      real*8  xi(6), xtol
      
      data xi/1.d0,1.d0,1.d0,0.d0,0.d0,0.d0/
      xtol = 1.d-12

c...  initialize material parameters 
      lam   = props(1)    ! lame constant
      mu    = props(2)    ! shear modulus
      xn0(1) = props(3)   ! n0[1]
      xn0(2) = props(4)   ! n0[2]
      xn0(3) = props(5)   ! n0[3]
      alpha = props(6)    ! growth rate
      tcr   = props(7)    ! critical stretch
      tmax  = props(8)    ! maximum growth
      gam   = props(9)    ! exponential growth rate

c...  calculate deformed elastic normal 
      norm = sqrt(xn0(1)*xn0(1) + xn0(2)*xn0(2) + xn0(3)*xn0(3))
      xn(1) = (dfgrd1(1,1)*xn0(1) + dfgrd1(1,2)*xn0(2) + dfgrd1(1,3)*xn0(3))/norm
      xn(2) = (dfgrd1(2,1)*xn0(1) + dfgrd1(2,2)*xn0(2) + dfgrd1(2,3)*xn0(3))/norm
      xn(3) = (dfgrd1(3,1)*xn0(1) + dfgrd1(3,2)*xn0(2) + dfgrd1(3,3)*xn0(3))/norm

      statev(4) = xn(1)
      statev(5) = xn(2)
      statev(6) = xn(3)

c...  calculate determinant of deformation gradient
      detf = +dfgrd1(1,1)*(dfgrd1(2,2)*dfgrd1(3,3)-dfgrd1(2,3)*dfgrd1(3,2))
     #       -dfgrd1(1,2)*(dfgrd1(2,1)*dfgrd1(3,3)-dfgrd1(2,3)*dfgrd1(3,1))
     #       +dfgrd1(1,3)*(dfgrd1(2,1)*dfgrd1(3,2)-dfgrd1(2,2)*dfgrd1(3,1))

c...  calculate area stretch
      finv(1,1) = (+dfgrd1(2,2)*dfgrd1(3,3) - dfgrd1(2,3)*dfgrd1(3,2))/detf
      finv(1,2) = (-dfgrd1(1,2)*dfgrd1(3,3) + dfgrd1(1,3)*dfgrd1(3,2))/detf
      finv(1,3) = (+dfgrd1(1,2)*dfgrd1(2,3) - dfgrd1(1,3)*dfgrd1(2,2))/detf
      finv(2,1) = (-dfgrd1(2,1)*dfgrd1(3,3) + dfgrd1(2,3)*dfgrd1(3,1))/detf
      finv(2,2) = (+dfgrd1(1,1)*dfgrd1(3,3) - dfgrd1(1,3)*dfgrd1(3,1))/detf
      finv(2,3) = (-dfgrd1(1,1)*dfgrd1(2,3) + dfgrd1(1,3)*dfgrd1(2,1))/detf
      finv(3,1) = (+dfgrd1(2,1)*dfgrd1(3,2) - dfgrd1(2,2)*dfgrd1(3,1))/detf
      finv(3,2) = (-dfgrd1(1,1)*dfgrd1(3,2) + dfgrd1(1,2)*dfgrd1(3,1))/detf
      finv(3,3) = (+dfgrd1(1,1)*dfgrd1(2,2) - dfgrd1(1,2)*dfgrd1(2,1))/detf

      ftn0(1) = finv(1,1)*xn0(1)+finv(2,1)*xn0(2)+finv(3,1)*xn0(3)
      ftn0(2) = finv(1,2)*xn0(1)+finv(2,2)*xn0(2)+finv(3,2)*xn0(3)
      ftn0(3) = finv(1,3)*xn0(1)+finv(2,3)*xn0(2)+finv(3,3)*xn0(3)
      
      the = detf*sqrt(ftn0(1)*ftn0(1) + ftn0(2)*ftn0(2) + ftn0(3)*ftn0(3))

c...  obtain state variable history
      theg_n = statev(1)
      theg   = theg_n

c...  local newton iteration 
c...  ------------------------------------------------------------------
      nitl = 0

      phig = (the/theg) - tcr

      if (phig.gt.0) then
 200    continue      
        nitl = nitl + 1

        kg  =    alpha*((tmax-theg)/(tmax-1.d0))**gam
        dkg =   -gam*kg/(tmax-theg)
        phig  =  the/theg-tcr
        dphig = -the/theg**2.d0

        res  = theg - theg_n - kg*phig*dtime
        dres = 1.d0 - (kg*dphig + dkg*phig)*dtime

        theg = theg - res / dres

        if ((nitl.lt.20).and.(dabs(res).gt.xtol)) go to 200
        if (nitl.eq.20) print *, 'no local convergence! |r|=',dabs(res)
        
        fac = kg*dtime/detf/dres/theg ! coefficient for growth tangent

      else
        theg = theg_n
        fac = 0.d0  
      endif
c...  ------------------------------------------------------------------
c...  end local newton iteration 

c...  update state variables 
      statev(1) = theg
      statev(2) = the/theg
      statev(3) = the

c...  calculate elastic tensor Fe = 1\sqrt(theg)*F + (1-1/sqrt(theg))*n \otimes n_0
      fe(1,1) = dfgrd1(1,1)/sqrt(theg) + (1.d0 - 1.d0/sqrt(theg))*xn(1)*xn0(1)
      fe(1,2) = dfgrd1(1,2)/sqrt(theg) + (1.d0 - 1.d0/sqrt(theg))*xn(1)*xn0(2)
      fe(1,3) = dfgrd1(1,3)/sqrt(theg) + (1.d0 - 1.d0/sqrt(theg))*xn(1)*xn0(3)
      fe(2,1) = dfgrd1(2,1)/sqrt(theg) + (1.d0 - 1.d0/sqrt(theg))*xn(2)*xn0(1)
      fe(2,2) = dfgrd1(2,2)/sqrt(theg) + (1.d0 - 1.d0/sqrt(theg))*xn(2)*xn0(2)
      fe(2,3) = dfgrd1(2,3)/sqrt(theg) + (1.d0 - 1.d0/sqrt(theg))*xn(2)*xn0(3)
      fe(3,1) = dfgrd1(3,1)/sqrt(theg) + (1.d0 - 1.d0/sqrt(theg))*xn(3)*xn0(1)
      fe(3,2) = dfgrd1(3,2)/sqrt(theg) + (1.d0 - 1.d0/sqrt(theg))*xn(3)*xn0(2)
      fe(3,3) = dfgrd1(3,3)/sqrt(theg) + (1.d0 - 1.d0/sqrt(theg))*xn(3)*xn0(3)

      detfe = +fe(1,1) * (fe(2,2)*fe(3,3)-fe(2,3)*fe(3,2))
     #        -fe(1,2) * (fe(2,1)*fe(3,3)-fe(2,3)*fe(3,1))
     #        +fe(1,3) * (fe(2,1)*fe(3,2)-fe(2,2)*fe(3,1))
     
      lnJe = dlog(detfe)

c...  calculate left cauchy-green deformation tensor be = fe * fe^t
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
      ddsdde(1,1) = (lam + 2.d0*(lam*lnJe - mu))/detfe + 2.d0*stress(1)
      ddsdde(2,2) = (lam + 2.d0*(lam*lnJe - mu))/detfe + 2.d0*stress(2)
      ddsdde(3,3) = (lam + 2.d0*(lam*lnJe - mu))/detfe + 2.d0*stress(3)
      ddsdde(1,2) = (lam)/detfe
      ddsdde(1,3) = (lam)/detfe
      ddsdde(2,3) = (lam)/detfe
      ddsdde(1,4) = stress(4)
      ddsdde(2,4) = stress(4)
      ddsdde(3,4) = 0.d0
      ddsdde(4,4) = (lam*lnJe - mu)/detfe + (stress(1) + stress(2))/2.d0
      
      if (ntens.eq.6) then
        ddsdde(1,5) = stress(5)
        ddsdde(2,5) = 0.d0
        ddsdde(3,5) = stress(5)
        ddsdde(1,6) = 0.d0
        ddsdde(2,6) = stress(6)
        ddsdde(3,6) = stress(6)
        ddsdde(5,5) = (lam*lnJe - mu)/detfe + (stress(1) + stress(3))/2.d0
        ddsdde(6,6) = (lam*lnJe - mu)/detfe + (stress(2) + stress(3))/2.d0
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
      cg_ij(1) = (lam*lnJe - lam - mu)*xi(1) + mu*xn(1)*xn(1)
      cg_ij(2) = (lam*lnJe - lam - mu)*xi(2) + mu*xn(2)*xn(2)
      cg_ij(3) = (lam*lnJe - lam - mu)*xi(3) + mu*xn(3)*xn(3)
      cg_ij(4) = (lam*lnJe - lam - mu)*xi(4) + mu*xn(1)*xn(2)
      cg_ij(5) = (lam*lnJe - lam - mu)*xi(5) + mu*xn(1)*xn(3)
      cg_ij(6) = (lam*lnJe - lam - mu)*xi(6) + mu*xn(2)*xn(3)

      cg_kl(1) = the*xi(1) - detf*detf/the*ftn0(1)*ftn0(1)
      cg_kl(2) = the*xi(2) - detf*detf/the*ftn0(2)*ftn0(2)
      cg_kl(3) = the*xi(3) - detf*detf/the*ftn0(3)*ftn0(3)
      cg_kl(4) = the*xi(4) - detf*detf/the*ftn0(1)*ftn0(2)
      cg_kl(5) = the*xi(5) - detf*detf/the*ftn0(1)*ftn0(3)
      cg_kl(6) = the*xi(6) - detf*detf/the*ftn0(2)*ftn0(3)

c...  compile tangent
      do i = 1,ntens
        do j = 1,ntens 
          ddsdde(i,j) = ddsdde(i,j) + fac*cg_ij(i)*cg_kl(j)
        enddo
      enddo

c...  calculate strain energy
      sse = (lam*lnJe**2.d0 + mu*(be(1)+be(2)+be(3) - 3.d0 - 2.d0*lnJe))/2.d0

      return
      end
