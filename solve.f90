subroutine solve

!use pks 
use system
use results
use solver
use const

implicit none
real*8 tolerancia,criterio,check_Ka_alpha,check_ka_beta,checkb_Kai,KK0check,KKaAcheckplus,kkaBcheckmin
real*8 x1(7)
real*8 x1g(7)
real*8 checkresults
integer ier, i,newt, j, ii, jj
integer flag
integer ngrid,iii

criterio=1E-6!criterio pra la norma
tolerancia=1E-6!criterio pra la diferencia de concentraciones relativa

x1(1)=log(xmBalphainitial) ! B in alpha
x1g(1)=x1(1)  

x1(2) = log(xmNaalphainitial)      ! initial guess for xmNaalpha
x1g(2) = x1(2)

x1(3)=log(xmClalphainitial) ! Cl in alpha
x1g(3)=x1(3)

x1(4)=log(xmAbetainitial)     !xAbeta  inicial 
x1g(4)=x1(4) 

x1(5)=log(xmBbetainitial)     !xBbeta  inicial
x1g(5)=x1(5)

x1(6)=log(xmNabetainitial)     !xNabeta  inicial
x1g(6)=x1(6)

x1(7)=log(xmClbetainitial)     !xClbeta  inicial
x1g(7)=x1(7)

!x1(8) =log(xmHplusalpha)!log(xmHplusbetainitial)      ! initial guess for xmNaalpha
!x1g(8) = x1(8)

!x1(9) =log(xmOHminalpha)! log(xmOHminbetainitial)      ! initial guess for xmNaalpha
!x1g(9) = x1(9)

!x1(10) = log(fA_HB_alphainitial)      ! initial guess for xmNaalpha
!x1g(10) = x1(10)

!print*,'in:alpha',ratioalpha,exp(x1(1)),exp(x1(2)),exp(x1(3))!,exp(x1(10))
!print*,'in:beta',exp(x1(4)),exp(x1(5)),exp(x1(6)),exp(x1(7))!,exp(x1(8)),exp(x1(9))

xSolventalpha=1-xmAalpha*Ma*vp*vsol-exp(x1(1))*Mb*vp*vsol-exp(x1(2))*vneg*vsol-exp(x1(3))*vneg*vsol
xSolventalpha=xSolventalpha/(1+expmuHplus+expmuOhmin)
xmSolventalpha=xSolventalpha/vs

xSolventbeta=1-exp(x1(4))*Ma*vp*vsol-exp(x1(5))*Mb*vp*vsol-exp(x1(6))*vneg*vsol-exp(x1(7))*vneg*vsol
xSolventbeta=xSolventbeta/(1+expmuHplus+expmuOhmin)
xmSolventbeta=xSolventbeta/vs

fhb_A_alpha=0
fhb_A_beta=0
!print*,'Na_alfa', 'n_tot_alfa', 'EO/Na alfa', 'Na_beta-Na_alpha', 'n_tot_beta', 'EO/Na beta'
call call_kinsol(x1, x1g, ier)

print*,'out:alpha',ratioalpha,exp(x1(1)),exp(x1(2)),exp(x1(3))!,exp(x1(10))
print*,'out:beta',exp(x1(4)),exp(x1(5)),exp(x1(6)),exp(x1(7))!,exp(x1(8)),exp(x1(9))


!checkresults=0.
!checkresults=abs((x1(1)- x1(2)/(x1(1)+x1(2)))) ! Naalpha - Nabeta
!checkresults=checkresults+abs((n_totalpha- x1(3))/(n_totalpha + x1(3))) ! n_totalpha - n_totbeta
!checkresults=checkresults+abs((log(ratioEOAalpha) - x1(4))/(log(ratioEOAalpha) - x1(4))) ! ratio alfa - ratio beta

!print*, norma
!print*, checkresults
!stop


!if ((norma.lt.criterio).and.( checkresults.gt.tolerancia)) then ! encuentra solucion
if (norma.lt.criterio) then ! encuentra solucion
!print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!yes'
print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!yes!!!!!!!!!!!!!!!!!!!!!!!!!!'
print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!yes!!!!!!!!!!!!!!!!!!!!!!!!!!'
!print*,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!yes!!!!!!!!!!!!!!!!!!!!!!!!!!'
!stop
write(8000,*) exp(x1(1)),exp(x1(2)),exp(x1(3)), exp(x1(4)), exp(x1(5)),exp(x1(6)),exp(x1(7)),xmAalpha ,pKD,ratioalpha
    print*,'Grid Point OK',yes
    if(justone.eq.1)call endall


! Found a solution, save in arrays

! Solutions for next iteration
      
      xmBalphainitial=exp(x1(1))! Ratio inicial alpha 1:1
      xmNaalphainitial=exp(x1(2))
      xmClalphainitial=exp(x1(3))     !xNabeta  inicial

      xmAbetainitial=exp(x1(4)) ! for next iteration
      xmBbetainitial=exp(x1(5)) !xratioeobeta inicial 
      xmNabetainitial=exp(x1(6))
      xmClbetainitial=exp(x1(7))
! save arrays

      xmBalpha = exp(x1(1))
      xmNaalpha= exp(x1(2))
      xmClalpha= exp(x1(3))  

      xmAbeta = exp(x1(4))
      xmBbeta = exp(x1(5)) 
      xmNabeta= exp(x1(6))
      xmClbeta= exp(x1(7))    

      yes=yes+1 ! counter of sucessful solutions

      write(9999,*)yes, xmsolventalpha,xmsolventbeta,xmClalpha, xmClbeta
      write(9998,*)yes,faspol_a_alpha, fAspol_A_alpha, fasio_B_alpha,fasio_A_alpha
      write(9997,*)yes,faspol_B_beta, fAspol_A_beta, fasio_B_beta,fasio_A_beta
      write(9996,*)yes, xmNaalpha,xmNabeta, xmNaalpha-xmNabeta
      write(9995,*)yes, packconst, neutralconst

      write(9994,*)yes, fch_a_alpha, funas_A_alpha, fch_B_alpha,funas_b_alpha
      write(9993,*)yes, fch_a_beta, funas_A_beta, fch_B_beta,funas_b_beta


      arrayaddedNaCl(yes)=xmaddedNaCl/Na*1.d24
  

      arraymNa(1,yes)=xmNaalpha/Na*1.d24
      arraymNa(2,yes)=xmNabeta/Na*1.d24

      arraymCl(1,yes)=xmClalpha/Na*1.d24
      arraymCl(2,yes)=xmClbeta/Na*1.d24

      arraymcsal(1,yes)=arraymNa(1,yes)+arraymCl(1,yes) ! salt can be calculated from Cl- only
      arraymcsal(2,yes)=arraymNa(2,yes)+arraymCl(2,yes)

      arraymA(1,yes)=MA*xmAalpha/Na*1.e24 ! in units of molar monomers 
      arraymA(2,yes)=MA*xmAbeta/Na*1.e24

      arraymB(1,yes)=MB*xmBalpha/Na*1.e24  ! in units of molar monomers
      arraymB(2,yes)=MB*xmBbeta/Na*1.e24

      arraympoltot(1,yes)=arraymA(1,yes)+arraymB(1,yes)
      arraympoltot(2,yes)=arraymA(2,yes)+arraymB(2,yes)

      arrayratioBA(1,yes)=arraymB(1,yes)/arraymA(1,yes)
      arrayratioBA(2,yes)=arraymB(2,yes)/arraymA(2,yes)

endif ! found solution
return
end subroutine

