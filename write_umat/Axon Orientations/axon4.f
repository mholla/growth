c...  define fiber directions - concave up
      if (kinc.eq.1) then
        if (coords(1).le.LL/two) then
          cc(1) = zero
          cc(2) = zero
          cc(3) = zero     
        else
          cc(1) = LL
          cc(2) = zero
          cc(3) = zero
        end if 
        
        rr(1) = coords(1) - cc(1)
        rr(2) = coords(2) - cc(2)
        norm = sqrt(rr(1)**two + rr(2)**two)
        
        if ((rr(1).eq.zero).and.(rr(2).eq.zero)) then
          xn(1) = one
          xn(2) = zero
          xn(3) = zero
        else
          xn(1) =  rr(2)/norm
          xn(2) = -rr(1)/norm
          xn(3) =  zero
        endif     

        statev(5) = xn(1)
        statev(6) = xn(2)
        statev(7) = xn(3)

      else
        xn(1) = statev(5)
        xn(2) = statev(6)
        xn(3) = statev(7)
      endif
      