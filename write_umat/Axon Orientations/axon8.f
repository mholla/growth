c...  define fiber directions - randomly
      if (kinc.eq.1) then
        call random_number(a)
        call random_number(b)
        a = 2*a - 1
        b = 2*b - 1
        norm = sqrt(a**2 + b**2)
      
        xn(1) =  a/norm
        xn(2) =  b/norm
        xn(3) = zero

        statev(5) = xn(1)
        statev(6) = xn(2)
        statev(7) = xn(3)
      
      else
        xn(1) = statev(5)
        xn(2) = statev(6)
        xn(3) = statev(7)
      endif
      