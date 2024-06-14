
subroutine electroneutroalpha(electroa)
use const
use system
use results
implicit none
real*8 electroa
!!
electroa= -fch_A_alpha*Ma*xmAalpha +fch_B_alpha*Mb*xmBalpha+xmNaalpha-xmClalpha
electroa=electroa +xmHplusalpha - xmOHminalpha

end  subroutine 
