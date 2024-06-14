
subroutine muA(potquimA)
use const
use system
use results
implicit none
real*8 potquimA
!!
potquimA= log(xmAalpha*vsol)-log(xmAbeta*vsol)- chi*MA*(MA*(xmAalpha-xmAbeta)&
+MB*(xmBalpha-xmBbeta))+MA*(log(fch_A_alpha)-log(fch_A_beta))&
-packconst*Ma*vpol*vsol+Ma*neutralconst
!
!print*,'A1',log(xmAalpha*vsol)-log(xmAbeta*vsol)
!print*,'A2',MA*(log(fA_unas_alpha)-log(fA_unas_beta))
!print*,'A3',-packconst*Ma*vpol*vsol
!print*,'A4',+Ma*neutralconst

end  subroutine 
