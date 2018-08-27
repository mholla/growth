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

      if      (cmname(1:3).eq.'COR') then
        call umat_cortex(stress,statev,ddsdde,sse,time,dtime,props,dfgrd1,ntens,ndi,nstatv,nprops,noel,npt,kstep,kinc)
      else if (cmname(1:3).eq.'SUB') then
        call umat_subcortex(stress,statev,ddsdde,sse,dtime,props,dfgrd1,ntens,ndi,nstatv,nprops,noel,npt,kstep,kinc,coords)
      else 
        print *, "Please indicate a valid material name"
      end if
     
      return
      end
      
