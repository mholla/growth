c...  ------------------------------------------------------------------
      subroutine sdvini(statev,coords,nstatv,ncrds,noel,npt,layer,kspt)
c...  ------------------------------------------------------------------
      include 'aba_param.inc'

      dimension statev(nstatv)

      statev(1)=1.0d0   !the_g
      statev(2)=1.0d0   !the_e
      statev(3)=1.0d0   !the
      statev(4)=1.0d0   !stress_fiber
      statev(5)=1.0d0   !f0_1
      statev(6)=1.0d0   !f0_2
      statev(7)=1.0d0   !f0_3
      statev(8)=1.0d0   !f_1
      statev(9)=1.0d0   !f_2
      statev(10)=1.0d0  !f_3
     
      return
      end

