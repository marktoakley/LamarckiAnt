      subroutine num_to_char(count,ccount)

      use amhglobals,  only:SO
      implicit none
      integer count,i,j,k
      character ccount*1,c(1)

      ccount=' '

      if (iabs(count).gt.99999999) then
        write(SO,*) 'number too large to convert',count
        stop
      endif

      j=0
      do while(count.ne.0)
        j=j+1
        i=mod(count,10)
        c(j)=char(i+48)
        count=(count-i)/10
      enddo

      do k=1,j
        ccount(j-k+1:j-k+1)=c(k)
      enddo

      return
      end

