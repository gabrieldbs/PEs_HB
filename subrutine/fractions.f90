subroutine fractions(vectf,vectfrac)
use system
use const
use results

implicit none
integer i,ITERAHB,NITERAHB
real*16 vectf(6),vectfrac(9)
real*16 xmphiA,xmphiB
real*16 xmphiNa,xmphiCl!,xmphiHplus,xmphiOHmin
real*16 aa, bb, cc,auxA,auxB,AuxC,discriminant,quadPlus,quadMinus
real*16 betaA,betaB,alphaA,alphaB,deltaA,deltaB,Omega,gammaA
real*16 fAS_B,fas_A,fasio_A,fasio_B,fC_A,fC_B,fhb_A
real*16 fnc_A,fnc_B,neutroconst_beta, neutroconst_alpha
real*16 K0Dcheck,K0BClcheck,K0ANacheck,K0Acheckplus,k0bcheck,k0hbcheck
real*16 phase,xsolvent,xmHplus,xmOHmin,xmphiSolv
real*16 auxas,kaso,auxbc,eta

real*16 j11,j12,j22,j21,det,dx,dy
real*16 a, b, c, d, x, y,x_new,y_new,f1,f2
real*16 x0, y0, tol
integer max_iter,iter
  ! Valores iniciales
  x0 = 0.1
  y0 = 0.1
  tol = 1.0e-4
  max_iter = 100000

phase=vectf(6) ! 

xmphiA=vectf(1)
xmphiB=vectf(2)
xmphiNa=vectf(3) 
xmphiCl=vectf(4)
xmphiSolv=vectf(5)
!fHB_A=vectf(7)
gammaA=K0HB*xmphiA*Ma*vaa!
!fHB_A=0
xmHplus=xmphisolv*expmuHplus
xmOHmin=xmphisolv*expmuOHmin

!poner en las notas las constantes auxiliares  explicitamente 
betaA=K0A*xmphiSolv/xmHplus ! ec 17 con *
betaB=K0B*xmphiSolv/xmOHmin ! ec 20 con **
!alphaA=(vsal*(xmphisolv*vsol)**vsal)/(K0ANa*xmphiNa*vsal*vsol)! 29  diferencia con löas notas  termino de lasal en el paquint 
!alphaB=(vsal*(xmphisolv*vsol)**vsal)/(K0BCl*xmphiCl*vsal*vsol)!32    chequeae qu combiene 
alphaA=(vsal*(1.0)**vsal)/(K0ANa*xmphiNa*vsal*vsol)! 29  diferencia con löas notas  termino de lasal en el paquint 
alphaB=(vsal*(1.0)**vsal)/(K0BCl*xmphiCl*vsal*vsol)!32    chequeae qu combiene 



deltaA = 1/(alphaA+alphaA/betaA+1)  !page 74  exbottles 
deltaB = 1/(alphaB+alphaB/betaB+1) ! page 74 exbottle combinar manuscritpoa

eta=xmphib*mb/(xmphia*ma)
auxBC=1+(1./betaB)+(1./alphaB)
Kaso=vab*vsol*K0D*xmphia*Ma
auxAs=1. +betaA+betaA/alphaA

fnc_A=x0
fas_B=y0
conver=0
do iter = 1, max_iter

  f1=1-eta*fas_b-fnc_A*(1+betaA+betaA/alphaA)-gammaA*fnc_A*fnc_A
  f2=fas_b-Kaso*betaA*fnc_A*alphaB*deltaB*(1-fas_B)

  j11=-(1+betaA+betaA/alphaA)-gammaA*2*fnc_A
  j12=-eta
  j21=-Kaso*betaA*alphaB*deltaB*(1-fas_B)
  j22=1+Kaso*betaA*fnc_A*alphaB*deltaB
  
  det = j11 * j22 - j12 * j21

  dx = (f1 * j22 - f2 * j12) / det
  dy = (f2 * j11 - f1 * j21) / det


! Actualizar las variables
  x_new = fnc_A - dx
  y_new = fas_B - dy

      ! Verificar la convergencia
   if (abs(x_new - FNC_a) < tol .and. abs(y_new - FAS_b) < tol) then
    fnc_a = x_new
    FAS_b = y_new
   print* ,'convergio fractions',fnc_a , FAS_b
        !stop
    conver=1
    exit
   endif

  fnc_a = x_new
  fas_B=  y_new
 !i stop
enddo
if (conver==0)then
    print *, '!!!!!!!no convergio ', max_iter, "iteraciones."
endif

fas_a=fas_b*Mb*xmphib/(Ma*xmphia)!  fracition aso  POL A
fc_A=fnc_A*betaA
fhb_A=gammaA*fnc_A*fnc_A
fasio_A=fc_A/alphaA

fasio_B=(1-fas_B)*deltaB! fraction as ion pol b
fC_B = fasio_B*alphaB !ec 32  fraction charged pol 
fnc_B=fc_B/betaB


!print*,'fhb',fHB_A

vectfrac(1)= fas_B
vectfrac(2)= fas_A
vectfrac(3)= fasio_B
vectfrac(4)= fasio_A
vectfrac(5)= fnc_B
vectfrac(6)= fnc_A
vectfrac(7)= fC_B 
vectfrac(8)= fC_A
vectfrac(9)= fHB_A

!print*,'fhb , gamma',fHB_A,gammaA
!print*,'fas_A,fas_b',fas_A,fas_B
!print*,'fasio_A,fasio_b',fasio_A,fasio_B
!print*,'fnc_A,fnc_b',fnc_A,fnc_B
!print*,'fc_A,fc_b',fc_A,fc_B
!print*,'fhb_A',fHb_A

!print*,'TESTEOS DE CONSTANTES'
       K0Acheckplus= -log10( (xmHplus/xmphisolv)*(&             !!
       fc_A )/fNC_A*xsolbulk*1.0d24/(Na*vsol))- pKaA              !! esto era para chequear pkaA
       K0Bcheck=(xmOHmin/xmphiSolv)*(fC_B/fnC_b) - K0B !

     K0ANacheck=((1)**vsal)*fasio_A/(fC_a)/(xmphiNA*vsol) - K0ANa
     K0BClcheck=((1)**vsal)/(xmphiCl*vsol)/((fC_B)/fasio_B) - K0BCl
!     K0ANacheck=((xmphisolv*vsol)**vsal)*fasio_A/(fC_a)/(xmphiNA*vsol) - K0ANa
!     K0BClcheck=((xmphiSolv*vsol)**vsal)/(xmphiCl*vsol)/((fC_B)/fasio_B) - K0BCl
       
       K0Dcheck=fAS_b/(fC_A*fC_B*xmphia*Ma*vsol*vab)! - k0D!!chekeo asoc
       K0HBcheck=fHB_a/(fnC_A*fnc_A*xmphiA*MA*vAA) - K0HB!chekeo asoc
testconst_as= K0Dcheck-K0d
testconst_hb= K0HBcheck
if (phase==1) then
conver_b=conver
testconst_as_B= K0Dcheck-K0d
testconst_hb_b= K0HBcheck
endif
! print*,'testeoK, K0Acheckplus,k0Bcheck ',K0Acheckplus,k0Bcheck
! print*,'testeoK,K0ANacheck,k0BClchec',K0ANacheck,k0BClcheck
! test kod no esta bien
! print*,'teste,  K0HBcheck' ,  K0HBcheck!
! print*,'testeo,K0Dcheck' , K0Dcheck- K0d
!print*,'phase',phase
!
end subroutine

