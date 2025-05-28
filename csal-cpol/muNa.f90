
subroutine muNa(potquimNa)
use const
use system
use results
implicit none
real*8 potquimNa

potquimNa= log(xmNabeta*vsol) +packconst*vsol*vpos +neutralconst -log(xmNaalpha*vsol)

end  subroutine 
