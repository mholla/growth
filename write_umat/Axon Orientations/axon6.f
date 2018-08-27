c...  define fiber directions - tangentially
      if (kinc.eq.1) then
        xn(1) = one
        xn(2) = zero
        xn(3) = zero  

        statev(5) = xn(1)
        statev(6) = xn(2)
        statev(7) = xn(3)
      
      else
        xn(1) = statev(5)
        xn(2) = statev(6)
        xn(3) = statev(7)
      endif
      