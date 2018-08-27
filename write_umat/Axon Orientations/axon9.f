c...  define fiber directions - tangentially at the cortex, radially at the core, with a smooth transition
      if (kinc.eq.1) then
        if (coords(1).le.LL/two) then
          theta = -pi/two/HH*coords(2)
        else
          theta =  pi/two/HH*coords(2)
        end if 

      xn(1) = one*cos(theta)
      xn(2) = one*sin(theta)
      xn(3) = zero

      statev(5) = xn(1)
      statev(6) = xn(2)
      statev(7) = xn(3)
      
      else
        xn(1) = statev(5)
        xn(2) = statev(6)
        xn(3) = statev(7)
      endif
      