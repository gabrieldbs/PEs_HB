subroutine fractions(vectf,vectfrac)
use system
use const
use results
implicit none
integer i,ITERAHB,NITERAHB
real*16 vectf(8),vectfrac(9)
real*16 xmphiA,xmphiB
real*16 xmphiNa,xmphiCl!,xmphiHplus,xmphiOHmin
real*16 aa, bb, cc,auxA,auxB,AuxC,discriminant,quadPlus,quadMinus
real*16 betaA,betaB,alphaA,alphaB,deltaA,deltaB,Omega,gammaA
real*16 fAS_B,fas_A,fasio_A,fasio_B,fC_A,fC_B,fhb_A
real*16 fnc_A,fnc_B,neutroconst_beta, neutroconst_alpha
real*16 K0Dcheck,K0BClcheck,K0ANacheck,K0Acheckplus,k0bcheck,k0hbcheck
real*8 phase

NITERAHB =10
phase=vectf(8) ! identifica fase lo necesito para neutroconst

!DO ITERAHB=0, NITERAHb
! print*,first_it
!if (first_it<2)then ! en la pirmera iteration entra dos veces una para alpha y una para vbta
   fhb_A=0
!  neutralconst_al =0! 
!  neutralconst_be=0
   gammaA= 0!
!   first_it=first_it+1
!endif

!if(first_it>=2 .and. phase==1)then
!  print*,'alpha'
!  first_it=first_it+1
!  neutralconst=1/2*log(K0HB)*log(fch_A_alpha*fch_A_alpha*xmAalpha*Ma*vaa*vsol/fhb_A_alpha) !chequear
!  gammaA= K0HB*xmphiA*Ma*vaa*exp(-2*neutralconst)!neutroconst_alpha) !33
!elseif(first_it>=2 .and. phase==2) then
!  print*,'beta' 
!  first_it=first_it+1
!  neutralconst=1/2*log(K0HB)*log(fch_A_beta*fch_A_beta*xmAbeta*Ma*vaa*vsol/fhb_A_beta) !chequear
!  gammaA= K0HB*xmphiA*Ma*vaa*exp(-2*neutralconst)!neutroconst_beta) !33
!endif

xmphiA=vectf(1)
xmphiB=vectf(2)
xmphiNa=vectf(3) 
xmphiCl=vectf(4)
xmphiSolv=vectf(5)
xmphiHplus=vectf(6)
xmphiOHmin=vectf(7)

!poner en las notas las constantes auxiliares  explicitamente 
betaA=K0A*xmphiSolv/xmphiHplus ! ec 17 con *
betaB=K0B*xmphiSolv/xmphiOHmin ! ec 20 con **
alphaA=(vsal*(xmphisolv*vsol)**vsal)/(K0ANa*xmphiNa*vsal*vsol)! 29  diferencia con löas notas  termino de lasal en el paquint 
alphaB=(vsal*(xmphisolv*vsol)**vsal)/(K0BCl*xmphiCl*vsal*vsol)!32    chequeae qu combiene 

deltaA = 1/(alphaA+alphaA/betaA+1)  !page 74  exbottles 
deltaB = 1/(alphaB+alphaB/betaB+1) ! page 74 exbottle combinar manuscritpoa

Omega = 1/(alphaA*deltaA*alphaB*deltaB*vab*vsol*K0D*xmphiB*Mb) !! page  75 ex bottle c ombinar manuscrtos
!print*,'b', Omega ,alphaA,deltaA,alphaB,deltaB,vab,K0D,xmphiB,Mb
AuxC = (1-fHB_A)*Ma*xmphia/(Mb*xmphib) !  eta   ! ec 26 docs
AuxB = - (1+Omega+(1-fHB_A)*(xmphiA/xmphiB))*(xmphiB/xmphiA)
AuxA =  1.
discriminant = sqrt(auxB**2 - 4.0*auxA/auxC)
quadPlus=(- auxB + discriminant)/(2.0*auxA)
quadMinus=(- auxB - discriminant)/(2.0*auxA)

    if((quadPlus.ge.0.0).and.(quadPlus.le.1.0))then
      fAS_A=quadPlus
      print*, 'Positive discriminant used in the fraction calculation!!!'
    !  stop
      !print*,'QUADPLUS',iR,iZ
   endif
   if((quadMinus.ge.0.0).and.(quadMinus.le.1.0))then
      fAS_a=quadMinus !fraction aso Pol B
      print*,'QUADMINUS',quadMinus,auxB ,discriminant,auxA
   endif
   if(((quadMinus.ge.0.0).and.(quadMinus.le.1.0)).and.((quadPlus.ge.0.0).and.(quadPlus.le.1.0)))then
      print*, 'Both quadratic solutions are between 0 and 1!!!'
    !  stop
   endif
fas_b=fas_a*Ma*xmphia/(Mb*xmphib)!  fracition aso  POL A

fasio_A= (1-fAS_a-fHB_A)*deltaA  ! fraction as ion pol A
fasio_B=(1-fas_B)*deltaB! fraction as ion pol b

fC_A = fasio_A*alphaA !ec 29 fraction charged POL-A
fC_B = fasio_B*alphaB !ec 32  fraction charged pol 
fnc_A=fc_A/betaA ! ec 17  *  fraction uncharged Pol-A
fnc_B=fc_B/betaB ! ec 20  ** fracion Uncharged POl- B

fHB_A=gammaA*fc_A*fc_A ! ec 33 no


vectfrac(1)= fas_B
vectfrac(2)= fas_A
vectfrac(3)= fasio_B
vectfrac(4)= fasio_A
vectfrac(5)= fnc_B
vectfrac(6)= fnc_A
vectfrac(7)= fC_B 
vectfrac(8)= fC_A
vectfrac(9)= fHB_A

print*,'fhb , gamma',fHB_A,gammaA
print*,'fas_A,fas_b',fas_A,fas_B
print*,'fasio_A,fasio_b',fasio_A,fasio_B
print*,'fnc_A,fnc_b',fnc_A,fnc_B
print*,'fc_A,fc_b',fc_A,fc_B
print*,'fhb_A',fHb_A

!print*,'TESTEOS DE CONSTANTES'
       K0Acheckplus= -log10( (xmphiHplus/xmphisolv)*(&             !!
       fc_A )/fNC_A*xsolbulk*1.0d24/(Na*vsol))- pKaA              !! esto era para chequear pkaA
      
       K0Bcheck=(xmphiOHmin/xmphiSolv)*(fC_B/fnC_b) - K0B !
      
       K0ANacheck=((xmphisolv*vsol)**vsal)*fasio_A/(fC_a)/(xmphiNA*vsol) - K0ANa
       K0BClcheck=((xmphiSolv*vsol)**vsal)/(xmphiCl*vsol)/((fC_B)/fasio_B) - K0BCl
       
       K0Dcheck=fAS_a/(fC_A*fC_B*xmphib*Mb*vsol*vab)! - k0D!!chekeo asoc
       K0HBcheck=fHB_a/(fC_A*fC_A*xmphiA*MA*vsol*vAA) - K0HB!!chekeo asoc
   
 !      print*,'testeoK',K0Acheckplus,k0Bcheck
!       print*,'testeoK',K0ANacheck,k0BClcheck
      print*,'testeo' , K0Dcheck- K0d

!stop 
!enddo  !iter

end subroutine

