!potquim OH- alpha- beta ! i.e. dilute - coacervate phase, following eq.xx

subroutine muOH(potquimOH)
use const
use system
use results
implicit none
real*8 potquimOH
!!
potquimOH= log(xmOHminbeta*vsol)-log(xmOHminalpha*vsol)+packconst*vsol+neutralconst

end  subroutine 
