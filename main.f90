  !  ###############################################################################!     
  !   Molecular Theory Program 
  !     H  base del pss/pdma
  !###############################################################################
!use pks
use system
use const
use solver
use results
implicit none
integer i, j, iii,isal
real*8 logxmAalpha
real*8 logratioalpha
real*8 xposbulk,xnegbulk,xHplusbulk,xOHminbulk,Kw,xsalt
real*8 pKw,KAinput,Kbind,KBinput,KANainput,KBClinput,kte
real*8 kbcl,KANa
integer k, kk, kkk
integer kkkk, kkkkk,kkkkkk

print*, ' HB_PEsolution GIT Version: ', _VERSION

call readinput ! read input from file
call allocation
vs=vsol
vaa=1.
vab=1.
vpol=vpolcero/vsol
vneg=4./3.*pi*rsal**3/vsol !volume of anion in units of vsol
yes=0 ! es para  chequear si encuentra o no xalpha, xbeta


pKw=14.0
kW=10**(-pKw)
pOHbulk=pKw-pHbulk
cOHminbulk=10**(-pOHbulk)
cHplusbulk=10**(-pHbulk)

vp=vpol!
vsal=vsol/vsol! Test
vpos=vsal
vneg=vsal

do k=1,pKDp
do kk=1,pKANap
do kkk=1,pKBClp
do kkkk = 1,Map 
do kkkkk = 1,Mbp 
do kkkkkk = 1,pKHBp

pKHB = pKHBi + (pKHBf-pKHBi)*float(kkkkkk-1)/float(pKHBp)  !asoc Pol-A -Pol-A

pKD = pKDi + (pKDf-pKDi)*float(k-1)/float(pKDp)  !asoc Pol-A -- Pol-B
pKANa = pKANai + (pKANaf-pKANAi)*float(kk-1)/float(pKANap) !asoc Pol-A Na
pKBCl = pKBCli + (pKBCLf-pKBCli)*float(kkk-1)/float(pKBClp)! asoc Pol-B Cl
Ma = Mai + (Maf-Mai)*float(kkkk-1)/float(Map) ! #  chains Pol A
Mb = Mbi + (Mbf-Mbi)*float(kkkkk-1)/float(Mbp) !# Chains Pol B


KHB=10**(-pKHB)
KD=10**(-pKD)
KA=10**(-pKaA)
KB=10**(-pKaB)
KANa=10**(-pKANa)
KBCl=10**(-pKBCl)

!xposbulk=phi_sal  !NUEVO
!xnegbulk=phi_sal  !NUEVO
print*,' pKHB ,pKD,pKANa,pKBCl ,Ma,Mb,pkaA,pkaB'

print*, pKHB ,pKD,pKANa,pKBCl ,Ma,Mb,pKaA,pKaB

print*,'pkaA pKaB',pkAa,pkaB

print*,'vpol',vpol,'vsol',vsol,'vneg',vsal

print*,' pKD' ,pKD,' pkHB' ,pkHB

xHplusbulk = (cHplusbulk*Na/(1.0d24))*(vs)
xOHminbulk = (cOHminbulk*Na/(1.0d24))*(vs)

xmHplusalpha=xHplusbulk
xmOHminalpha=xOHminbulk
 ! do isal = 1, ncsal ! loop in csal

 ! csal = csalini + (csalfin-csalini)/float(ncsal-1)*float(i-1)


 ! xsalt=csal*vsal*vs*6.02/10.0

  ! if(pHbulk.le.7) then  ! pH<= 7
  !   xposbulk=xsalt
  !   xnegbulk= xsalt +(xHplusbulk -xOHminbulk)*vsal! NaCl+ HCl  
  !! else                  ! pH >7 
  !   xposbulk=xsalt +(xOHminbulk -xHplusbulk)*vsal ! NaCl+ NaOH   
  !   xnegbulk=xsalt
  ! endif

! Concentration of free anf paired Na+ and Cl- in bulk reference

!  rhoNa = xposbulk/vs
!  rhoCl = xnegbulk/vs

!  aa = 1.
!  bb = -1.*(rhoNa+rhoCl+1./Ksal/vs)
!  cc = rhoNa*rhoCl

!  rhoNaCl = ((-bb - sqrt(bb**2 - 4.*aa*cc))/(2.0*aa))

!  rhoNa = rhoNa - rhoNaCl
!  rhoCl = rhoCl - rhoNaCl

!  print*, rhoNa*rhoCl/rhoNaCl

!  xposbulk = rhoNa*vs
!  xnegbulk = rhoCl*vs
!  xNaClbulk = rhoNaCl*2.*vs

!  print*, '!!!!', xposbulk, xnegbulk, xNaClbulk
!  print*, '$$$$', (xposbulk/vs)*(xnegbulk/vs)/(xNaClbulk/2.0/vs), Ksal

!  if(xNaClbulk.lt.0.0)cycle ! if negative, go to next salt concentration
!!!!!

  xsolbulk=1.0 -xHplusbulk -xOHminbulk! - xnegbulk -xposbulk! -xNaClbulk

  K0ANa = (KANa)*(xsolbulk**vpos/vsol/(Na/1.0d24)) ! thermodynamic constants
  K0BCl = (KBCL)*(xsolbulk**vneg/vsol/(Na/1.0d24))
  K0D = (KD)*(Na/1.0d24)
  K0HB = (KHB)*(Na/1.0d24)
  K0A = (KA*vs/xsolbulk)*(Na/1.0d24)! intrinstic equilibruim constant 
  K0B = (Kw/KB*vs/xsolbulk)*(Na/1.0d24)
 
  xmsolventalpha=xsolbulk
  xmsolventbeta=xsolbulk
!  expmupos=xposbulk/xsolbulk**vsal
!  expmuneg=xnegbulk/xsolbulk**vsal
  expmuHplus=xHplusbulk/xsolbulk ! vHplus=vsol
  expmuOHmin=xOHminbulk/xsolbulk ! vOHminus=vsol
  
  do j=1, npasosratioalpha  ! loop over ratio Pol-A/PolB en  alpha
 
    logratioalpha = (logratioalphaf-logratioalphai)*float(j-1)/float(npasosratioalpha) + logratioalphai
    ratioalpha= 10**(logratioalpha) !10**  !Segunda variable que fijamos  xmpoltotalalpha


    do i = 1, npasosxmAalpha ! loop over Pol-A en alpha

      first_it= 0 ! para fHB pårimero 0  y despues itera
      logxmAalpha = logxmAalphai  + (logxmAalphaf-logxmAalphai) &
      /float(npasosxmAalpha)*float(i-1)  !Na
 
      xmAalpha = 10**(logxmAalpha)

      iter=0
      print*,'yes'
      call solve

     enddo ! j
  
  enddo ! i
 
!  enddo  !isal

enddo !k
enddo !kk
enddo !kkk
enddo !kkkk
enddo !kkkkk
enddo !kkkkkk


! SAVE RESULTS TO FILE

open (unit=3,file='csal_poltot_mol_alpha.txt',status='replace')

do iii=1,yes
   write (3,*) arraympoltot(1,iii), arraymcsal(1,iii)
end do

open (unit=4,file='csal_poltot_mol_beta.txt',status='replace')

do iii=1,yes
   write (4,*) arraympoltot(2,iii), arraymcsal(2,iii)
end do

open (unit=40,file='cpoltot_ratioBA_mol_alpha.txt',status='replace')

do iii=1,yes
   write (40,*) arrayratioBA(1,iii), arraympoltot(1,iii)
end do

open (unit=30,file='cpoltot_ratioBA_mol_beta.txt',status='replace')

do iii=1,yes
   write (30,*) arrayratioBA(2,iii), arraympoltot(2,iii)
end do

open (unit=400,file='csal_ratioBA_mol_alpha.txt',status='replace')

do iii=1,yes
   write (400,*) arrayratioBA(1,iii), arraymcsal(1,iii)
end do

open (unit=300,file='csal_ratioBA_mol_beta.txt',status='replace')

do iii=1,yes
   write (300,*) arrayratioBA(2,iii), arraymcsal(2,iii)
end do

open (unit=600,file='polA_polB_alpha.txt',status='replace')

do iii=1,yes
   write (600,*) arraymA(1,iii), arraymB(1,iii)
end do
 

open (unit=500,file='polA_polB_beta.txt',status='replace')

do iii=1,yes
   write (500,*) arraymA(2,iii), arraymB(2,iii)
end do

open (unit=600,file='polA_addedNa_alpha.txt',status='replace')

do iii=1,yes
   write (600,*) arraymA(1,iii), arrayaddedNaCl(iii)
end do



call endall     ! clean up and terminate
end 



subroutine endall
 stop
end subroutine
