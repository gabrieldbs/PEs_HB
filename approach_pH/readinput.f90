subroutine readinput

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!  This routine reads variables from fort.8
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!use pks
use system
use solver
!use kai

implicit none
integer i
real*8 xmNaalphaini
character*30 basura1
character*30 basura2
character*30 basura3
character*30 basura4
character*30 basura5
character*30 basura6
character*30 basura7
character*30 basura8
character*30 basura9
character*30 basura10
character*30 basura11
character*30 basura12
character*30 basura13
character*30 basura14
character*30 basura15

! read starts here, not that read is performed sequentially! 

read(8,*)  basura1
read(8,*)  Mai, Maf, Map    ! Ma  polA

read(8, *)  basura2
read(8, *)  Mbi, Mbf, Mbp  !   MEO  polEO

read(8, *) basura3! vp
read(8, *)  vpolcero !

read(8, *)  basura4! rsal
read(8, *)  rsal !

read(8, *)  basura5
read(8, *)  logratioalphai, logratioalphaf, npasosratioalpha  !  scan total monomer density in molar

!read(8, *) basura5
!read(8, *) csalini,csalfin,ncsal  ! salt concentration in bulk (Molar)

read(8, *)  basura6
read(8, *)  logxmAalphai, logxmAalphaf, npasosxmAalpha   ! scan ratio monomer density
!
read(8, *)  basura7
read(8, *)  pKDi, pKDf, pKDp     ! polymer-polymer attraction strenght in kBT
!
read(8, *)  basura8
read(8, *) pKaA     ! polymer-polymer attraction strenght in kBT
!
read(8, *) basura9
read(8, *) pkaB     ! polymer-polymer attraction strenght in kBT
!
read(8, *)  basura10
read(8, *) pKANai, pKANaf, pKANap      ! polymer-polymer attraction strenght in kBT
!
read(8, *) basura11
read(8, *) pkBCli, pKBClf, pKBClp     ! polymer-polymer attraction strenght in kBT

read(8, *) basura12
read(8, *) pkHBi, pKHBf, pKHBp     ! polymer-polymer attraction strenght in k

read(8,*) basura12
read(8,*) pHbulk

read(8, *) basura13
read(8, *) chi  ! cutoff for3porr sv interaction in lattice sites


read(8,*) basura14
read(8,*) justone ! solves only one point

read(8, *) basura15
read(8, * ) xmBalphainitial,xmNaalphainitial,xmClalphainitial, &
            xmAbetainitial,xmBbetainitial, xmNabetainitial,xmClbetainitial ,  &
            xmHplusbetainitial,xmOHminbetainitial,fA_HB_alphainitial! read initial guess in the same order as output, ratio and xmtot in alpha are not solved for, so the result is stores in basura

end subroutine
