c...  define fiber directions - tangentially near the cortex, and radially at the core
      if (kinc.eq.1) then
        if (coords(2).ge.-HH/4.d0) then     ! tangential
          xn(1) = one
          xn(2) = zero
          xn(3) = zero
        else                                ! radial
          xn(1) = zero
          xn(2) = one
          xn(3) = zero
        end if 

      statev(5) = xn(1)
      statev(6) = xn(2)
      statev(7) = xn(3)
      
      else
        xn(1) = statev(5)
        xn(2) = statev(6)
        xn(3) = statev(7)
      endif
      