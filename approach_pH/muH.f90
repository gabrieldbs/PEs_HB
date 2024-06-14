
subroutine muH(potquimH)
use const
use system
use results
implicit none
real*8 potquimH
!!
potquimH= log(xmHplusbeta*vsol)-log(xmHplusalpha*vsol)+packconst*vsol+neutralconst

end  subroutine 
