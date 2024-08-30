subroutine fkfun(x,f,ier)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! User provided routine for kinsol
! x is the input vector
! f is the output vector, kinsol will change x in order to get f = 0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


use solver
use system
use const
use results
implicit none
 
real*8 xmsalt
real*8 Penality,testKa,testkeo,testkd
real*8 testneuta,testneutb,testpcka,testpckb
integer*4 ier
real*16 vectfalpha(6),vectfbeta(6),vectfrac(9)
real*8 x(7),f(7)
real*8 potA,potB,potNa,free_ener,potH,potOH,potquimH,potquimOH
real*8 muAalpha,muAbeta,muBalpha,muBbeta,fealpha,febeta
real*8 potquimA,elib,potquimB,potquimNa,muNaalpha,muNabeta
real*8 neutralalpha,neutralbeta,elecneualpha,elecneubeta
real*8 penalityA,penalityNa,penalityB
real*8 diffNaalpha
integer i


!xmAalpha
!xmHplusalpha=xHplusbulk
!xmOHalpha=xOHminbulk
xmBalpha =xmAalpha/ratioalpha  !b alpha demas para convergencia
xmNaalpha=exp(x(2)) !na alphga
xmClalpha= exp(x(3)) !cl alpha

xmAbeta=exp(x(4))        !xmA en beta
xmBbeta=exp(x(5))       !B beta
xmNabeta=exp(x(6))        !xmNa en beta
xmClbeta=exp(x(7))      !cl beta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

vectfalpha(1)=xmAalpha
vectfalpha(2)=xmBalpha
vectfalpha(3)=xmNaalpha
vectfalpha(4)=xmClalpha
vectfalpha(5)=xmSolventalpha
xmHplusalpha=xmSolventalpha*expmuhplus
xmOHminalpha=xmSolventalpha*expmuohmin
vectfalpha(6)=1
!vectfalpha(7)=fHB_A_alpha

vectfbeta(1)=xmAbeta
vectfbeta(2)=xmBbeta
vectfbeta(3)=xmNabeta
vectfbeta(4)=xmClbeta
vectfbeta(5)=xmSolventbeta
xmHplusbeta=xmSolventbeta*expmuhplus
xmOHminbeta=xmSolventbeta*expmuohmin
vectfbeta(6)=2
!vectfbeta(7)=fHB_A_beta


call fractions(vectfalpha,vectfrac)
!!!Aca  DEFINIRR BIEN LO VOY A NECESITAR
!vectfractions-> (1) fEO_aspol,(2)fA_aspol,(3) fEO_asion (4)fA_Asion,(5)fEO_unas,f(6)fA_unas


faspol_B_alpha=vectfrac(1)
faspol_A_alpha=vectfrac(2)

fasio_B_alpha=vectfrac(3)
fasio_A_alpha=vectfrac(4)

funas_B_alpha=vectfrac(5)
funas_A_alpha=vectfrac(6)

fch_B_alpha=vectfrac(7)
fch_A_alpha=vectfrac(8)

fHB_A_alpha=vectfrac(9)


call fractions(vectfbeta,vectfrac)

!vectdddfractions-> (1) fEO_aspol,(2)fA_aspol,(3) fEO_asion (4)fA_Asion,(5)fEO_unas,f(6)fA_unas

faspol_B_beta=vectfrac(1)
faspol_A_beta=vectfrac(2)

fasio_B_beta=vectfrac(3)
fasio_A_beta=vectfrac(4)

funas_B_BEta=vectfrac(5)
funas_A_beta=vectfrac(6)

fch_B_beta=vectfrac(7)
fch_A_beta=vectfrac(8)

fHB_A_beta=vectfrac(9)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  AUXILIARY CALC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!xmClalpha=-fA_unas_alpha*Ma*xmAalpha + fB_asion_alpha*Meo*xmEOalpha + xmNaalpha  ! xmCl en alpha
!xmClbeta=-fA_unas_beta*Ma*xmAbeta + fB_asion_beta*MB*xmEObeta + xmNabeta  ! xmCl en alpha


xSolventalpha=1. -Ma*vpol*vsol*xmAalpha -Mb*vpol*vsol*xmBalpha&
!-fasio_A_alpha*Ma*xmAalpha*vpos*vsol-fasio_B_alpha*Mb*xmBalpha*vneg*vsol&
-xmNaalpha*vpos*vsol -xmClalpha*vneg*vsol-(xmhplusalpha+xmohminalpha)*vsol

xSolventbeta=1. - Ma*vpol*vsol*xmAbeta  -Mb*vpol*vsol*xmBbeta&
!-fasio_A_beta*Ma*xmAbeta*vpos*vsol-fasio_B_beta*Mb*xmBbeta*vneg*vsol&
-xmNabeta*vpos*vsol -xmClbeta*vneg*vsol-(xmhplusbeta+xmohminbeta)*vsol

xmSolventalpha=xSolventalpha/vsol
xmSolventbeta =xSolventbeta/vsol

xmHplusalpha=xmSolventalpha*expmuhplus
xmOHminalpha=xmSolventalpha*expmuohmin
xmHplusbeta=xmSolventbeta*expmuhplus
xmOHminbeta=xmSolventbeta*expmuohmin

packconst=(1./vsol)*(log(xSolventalpha)-log(xSolventbeta) ) ! betapi 
neutralconst=log(xmClbeta*vsol)+vneg*vsol*packconst-log(xmClalpha*vsol) ! phi
!print*,'neutralconst', neutralconst,xmClbeta,xmClalpha,packconst

!print*,'Solven_alpha, Solven_beta,packconst,neutralconst',xSolventalpha,xSolventbeta,packconst,neutralconst
!!print*,' neutralconst',neutralconst,neutralconstalpha,neutralconstbeta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
!  CHEMICAL POTS AND FREE ENERGY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

potquimNa=0.
call muNa(potquimNa)
potNa=potquimNa

potquimA=0.
call muA(potquimA)
potA = potquimA

potquimB=0.
call muB(potquimB)
potB = potquimB

elib=0.
call fe(elib)
free_ener=elib

elecneualpha=0.
call electroneutroalpha(elecneualpha)
neutralalpha=elecneualpha

elecneubeta=0.
call electroneutrobeta(elecneubeta)
neutralbeta=elecneubeta

!potquimH=0.
!call muH(potquimH)
!potH=potquimH

!potquimOH=0.
!call muOH(potquimOH)
!potOH=potquimOH


xmaddedNaCl = xmNaalpha + fasio_A_alpha*MA*xmAalpha - xmAalpha*MA

Penality=abs(xmNabeta-xmNaalpha)/(xmNabeta*0.5+xmNaalpha*0.5)
Penality=Penality+abs(xmAalpha-xmAbeta)/(xmAalpha*0.5+xmAbeta*0.5)
Penality=Penality+abs(xmBalpha-xmBbeta)/(xmBalpha*0.5+xmBbeta*0.5)
Penality=Penality+abs(xmClalpha-xmClbeta)/(xmClalpha*0.5+xmClbeta*0.5)

!Penality=Penality+abs(ratioBAalpha-ratioEOAbeta)/(ratioeoaalpha*0.5+ratioeoabeta*0.5)

f(1)=-free_ener/Penality
f(2)=-potNa/Penality
f(3)=-potA/Penality
f(4)=-potB/Penality
f(5)=-neutralalpha/Penality
f(6)=-neutralbeta/Penality
f(7)=exp(x(1))-xmBalpha
!f(8)=exp(x(8))-xmHplusalpha!-potH/penality
!f(9)=exp(x(9))-xmOHminalpha !-potOH/penality
!f(10)=exp(x(10))-fhB_A_alpha
iter = iter + 1
norma = 0.0

do i = 1, 7
!print*,'f',i,f(i)
  norma = norma +(f(i))**2  
enddo

!if (norma<1.0E-002 )then
!      print*,'yes',xmaalpha,xmbalpha,xmNaalpha,xmClalpha
!      print*,'yes',xmabeta,xmbbeta,xmNabeta,xmClbeta
!stop
!endif

!testneuta=xmNaalpha -xmClalpha-Ma*xmAalpha*fA_unas_alpha+fEO_asion_alpha*Meo*xmEOalpha
!testneutb=xmNabeta -xmClbeta-Ma*xmAbeta*fA_unas_beta+fEO_asion_beta*Meo*xmEObeta

!testpcka=1.-Ma*xmAalpha*vpol*vsol-xmSolventalpha*vsol -Meo*xmEOalpha*vpol*vsol-xmNaalpha*vpos*vsol-xmClalpha*vneg*vsol
!testpcka=testpcka -vneg*vsol*(Ma*xmAalpha*(fA_asion_alpha+fA_aspol_alpha)+Meo*xmEOalpha*fEO_asion_alpha)
 
!testpckb=1.-Ma*xmAbeta*vpol*vsol-xmSolventbeta*vsol -Meo*xmEObeta*vpol*vsol-xmNabeta*vpos*vsol-xmClbeta*vneg*vsol
!testpckb=testpckb -vneg*vsol*(Ma*xmAbeta*(fA_asion_beta+fA_aspol_beta)+Meo*xmEObeta*fEO_asion_beta)
!print*,'testneutra',testneuta,testneutb,testpcka,testpckb

!if( iter==10) then
!stop
!endif
print*,norma,penality !if (Norma.lt.1E-5)then
!stop
ier = 0.0
return
end subroutine
