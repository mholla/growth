c...  define fiber directions - radially
      if (kinc.eq.1) then
        xn(1) = zero
        xn(2) = one
        xn(3) = zero  

        statev(5) = xn(1)
        statev(6) = xn(2)
        statev(7) = xn(3)
      
      else
        xn(1) = statev(5)
        xn(2) = statev(6)
        xn(3) = statev(7)
      endif
      