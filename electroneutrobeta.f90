
subroutine electroneutrobeta(electrob)
use const
use system
use results
implicit none
real*8 electrob
!!
electrob= -fch_A_beta*Ma*xmAbeta +fch_B_beta*Mb*xmBbeta+xmNabeta-xmClbeta
electrob= electrob +xmHplusbeta- xmOHminbeta

end  subroutine 
