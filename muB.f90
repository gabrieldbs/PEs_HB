
subroutine muB(potquimB)
use const
use system
use results
implicit none
real*8 potquimB
!!chequear eleccontstraini
potquimB=log(xmBalpha*vsol)-log(xmBbeta*vsol)-chi*MB*(Ma*(xmAalpha-xmAbeta)&
+MB*(xmBalpha-xmBbeta))+MB*(log(funas_B_alpha)-log(funas_B_beta))&
-packconst*MB*vpol*vsol-MB*neutralconst

!print*,'B1',log(xmBalpha*vsol)-log(xmBbeta*vsol)
!print*,'B2',MB*(log(fB_unas_alpha)-log(fB_unas_beta))
!print*,'B3',-packconst*MB*vpol*vsol
!print*,'B4',-MB*neutralconst




end subroutine 
