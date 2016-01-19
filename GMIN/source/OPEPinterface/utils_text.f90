subroutine convert_to_chain(init_number,chain, caractere)
  integer, intent(in) :: init_number
  character(len=30), intent(out) :: chain
  character(len=10) :: digits = '0123456789'
  character(len=1), intent(in) :: caractere

  integer :: i, decades, divider, remainder,number

  ! We first fill the chain with zeroes
  do i = 1, len(chain)
    chain(i:i) = caractere
  end do

  number = init_number
  if (number .le. 0 ) then
     chain(30:30) = "0"
     return
  endif
  
  decades = log10( 1.0d0 * number) + 1

  divider = 1
  do i=2, decades
    divider =  divider * 10 
  enddo

  do i = len(chain)-decades+1, len(chain)
    remainder = number / divider  + 1
    chain(i:i) =  digits(remainder:remainder)
    remainder = remainder -1
    number = number - remainder * divider
    divider = divider / 10
  enddo

end subroutine     

subroutine real_to_chain(init_number,width,decimal,chain)
  double precision, intent(in) :: init_number
  integer, intent(in) :: width, decimal

  character(len=30) :: textstring
  character(len=width) :: chain
  character(len=1) :: signe

  integer :: i, number, irest
  double precision :: rest

  if (init_number .lt. 0 ) then
     signe = "-"
     number = -1.0 * init_number
     rest  = -1.0*init_number - number
  else
     signe = " "
     number = init_number
     rest  = init_number - number
  endif


  do i = 1, decimal
     rest = rest*10
  end do
  irest = rest

  ! We first fill the chain with zeroes
  do i = 1, width
    chain(i:i) = " "
  end do

  ! Then we convert both chains separately
  call convert_to_chain(number,textstring," ")

  ! We decide whether or not we must add a minus sign
  do i=2, len(textstring) 
     if ( textstring(i:i) .ne. " " ) then
        textstring(i-1:i-1) = signe
        exit
     endif
  end do

  min = 1
  max = width - (decimal-1)
  chain(min:max) = textstring(31-(width-decimal-1):30)
  
  min = width-decimal
  max = min
  chain(min:max) = "."

  call convert_to_chain(irest,textstring,"0")
  min = width-decimal+1
  max = width
  chain(min:max) = textstring(31-decimal:30)

  return
end subroutine real_to_chain
