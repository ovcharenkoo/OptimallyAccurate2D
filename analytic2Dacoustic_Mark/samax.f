c
c
c  function to return the absolute maximum of an array
c
      function samax(n,x,inc)
      implicit none
      real*8 smx,samax
        
      dimension  x(n)
c
      smx = 0.0
          do 1100 i = 1,n,inc
          smx = max(abs(x(i)),smx)
1100      continue
      samax = smx
c
      return
      end
