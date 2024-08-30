module system
!real*16 delta   ! delta is the discretization lenght in z direction
!real*16 sigmaA
integer Ma,MB
real*8 rsal
integer ntot
real*8 vpolcero,vaa!real*8 xmNabetainitial, deltaxmNabetainitial, n_totbetainitial,ratioEOAbetainitial
real*8 xmBalphainitial,xmClalphainitial,xmAbetainitial,xmNabetainitial,xmBbetainitial,xmClbetainitial
real*8 vpol,vpos,vneg,vab
real*8 phimin, phimax
real*8 csal
real*8 chi
integer justone
integer yes
integer  first_it
real*8 pKA,pKB,pkD,pKHB,pkAna,pKBcl
real*8 pKAi,pKBi,pkDi,pKHBi
real*8 pKAf,pKBf,pkDf,pkHBf
real*8 pKANai,pKBCli
real*8 pKANaf,pKBClf
integer pKAp,pKBp,pkDp,pkHBp,pkANap,pkBclp
real*8 Mai, Maf, Mbi, Mbf
integer Map, Mbp
real*8 K0A, K0B,K0D,Ka,KB,Kd,K0HB,KHB
real*8 logratioalphai, logratioalphaf, logxmAalphai,logxmAalphaf
real*8 xmaddedNaCl
integer npasosxmAalpha, npasosratioalpha 
real*8 cNaplus,cClmin
real*8 xNaalpha,xmNaalpha,xmNabeta, xmNaalphatot
real*8 n_totalpha, n_totbeta
real*8 xmNaalphainitial
real*8 xClalpha,xClbeta,xmClalpha,xmClbeta
real*16 xmohminalpha,xmohminbeta,xmhplusalpha,xmhplusbeta
real*8 xsolventalpha,xsolventbeta
real*8 xmsolventalpha,xmsolventbeta
real*8 xAalpha,xAbeta,xBalpha,xBbeta
real*8 ratioalpha, xmAalpha,xmAbeta,xmBalpha,xmBbeta
real*8 ratioBAalpha,ratioBAbeta 
real*8 xmsalalpha,xmsalbeta
integer npasosratioEO
real*8 logratioEOAalphai, logratioEOAalphaf
real*8 vsal,vs,vp
real*16 pKaHplus,pKHbulk,cHplusbulk,pHbulk
real*16 pKaOHmin,pOHbulk,cOHminbulk
real*16 expmupos, expmuneg, expmuHplus, expmuOHmin
integer ncsal
real*16 K0ANa,K0BCl
real*16 pKaA
real*16 pkaB
real*8 xmphiOHmin,xmphiHplus,xsolbulk
real*8 fA_HB_alphainitial,xmhplusbetainitial,xmohminbetainitial
endmodule

module results
integer conteo, flaggg
real*8 faspol_A_alpha,fasio_A_alpha,funas_A_alpha
real*8 faspol_B_alpha,fasio_B_alpha,funas_B_alpha
real*8 faspol_A_beta,fasio_A_beta,funas_A_beta
real*8 faspol_B_beta,fasio_B_beta,funas_B_beta
real*8 testconst_as, testconst_hb,testconst_as_B, testconst_hb_b
real*8 fHB_A_alpha,fHB_A_beta
real*8 fch_A_alpha,fch_A_beta
real*8 fch_B_alpha,fch_B_beta
real*8 neutralconstalpha, neutralconstbeta
real*8 packconst,neutralconst
real*8 arraympoltot(2,100000)
real*8 arraymcsal(2,100000)
real*8 arraypoltot(2,100000)
real*8 arraycsal(2,100000)

real*8 arraymNa(2,100000)
real*8 arrayaddedNaCl(100000)
real*8 arraymCl(2,100000)
real*8 arraymA(2,100000)
real*8 arraymB(2,100000)
real*8 arrayratioBA(2,100000)
integer conver,conver_b
integer cont
endmodule

module const
real*8, parameter :: pi = 3.14159 ! pi 
real*8, parameter :: Na = 6.02d23 ! Avogadro's number
real*8, parameter :: vsol = 0.03  ! bjerrum lenght in water in nm
!real*16 constq
!real*16 pKw
endmodule


module solver
real*8 norma
integer iter
endmodule


module mkinsol
double precision, allocatable :: pp(:)
endmodule
