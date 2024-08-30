
long="50"
vpol="0.15"
pK0="-1.2 -1.3 -1.4 -1.4 -1.5 -1.6 -1.7 -1.8 -2"
pKAio="-1. -1.2 -1.5 -2"
pkHB="-1. -1.2 -1.5 -2 "
pH="7"
chi="0 0.1 0.15 0.2 0.3 0.5 0.7 1"

for j in $long ; do
echo $j

mkdir 29ag_$j
cd 29ag_$j 

for pHs in $pH ; do
echo $pHs

mkdir pH_$pHs
cd pH_$pHs

for vp in $vpol ; do
echo $vp

mkdir vpol_$vp
cd vpol_$vp

for pK in $pK0 ; do
echo $pK

mkdir pk_$pK
cd  pk_$pK

for pKio in $pKAio ; do
echo $pKio

mkdir pkion_$pKio
cd  pkion_$pKio

for pHB in $pkHB ; do
echo $pHB

mkdir pKHB_$pHB
cd  pKHB_$pHB

for cchi in $chi ; do
echo $cchi

mkdir chi_$cchi
cd  chi_$cchi

echo "
# Ma#
$j $j 1
#  Mb #
$j $j 1
#vpol
$vp
#rsal
0.27
# log10(ratio) min max npasos 3
0. 1 1 #1 40
#log10(xma) min max npasos
-6.49 1 300
#pk0#
$pK $pK 1 
#pkaA
4.5
#pKAB
9
#pkA_Aio
$pKio $pKio 1
#pkA_Bio
$pKio $pKio 1
#pkA_Bio
$pHB $pHB 1
#pH
$pHs
#chi#
$cchi 
#juiston0#
0
#xmAinitalpha(basura) xmNainitalpha(basura)  xmBinitalpha xmClinitalph xmAinitbeta xmNainitbeta  xmBinitbeta xmClinitbeta#
3.2359365692962813E-006   1.1269574282838625E-001 0.142123815043626389    1.0465340699255346E-002   1.2086845621600142E-002  1.1468290817675662E-001  0.2052032679637574 0.0001 0.0001 0.00001
$infile" > fort.8
~/Desktop/proyects/PEs_HB/subrutine/PEs-HB

cd ..

done
cd ..

done

cd ..

done
cd ..

done
cd ..

done

cd ..

done
cd ..

done

