#!/bin/bash
rm -r -f DATA_testn_part DATA_testn extract.dat extractpart.dat
mkdir DATA_testn_part
sed -e "s/XPARTX/1/" tests_intpartdistplan.par >tests_intpartstatntmp.par
../intpartstat.x  tests_intpartstatntmp.par >/dev/null || exit 1
mv DATA_testn_part DATA_testn

mkdir DATA_testn_part

for group in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 
do
echo "integration de 17 particules par groupe de" $group
rm -f DATA_testn_part/*
sed -e "s/XPARTX/$group/" tests_intpartdistplan.par >tests_intpartstatntmp.par
../intpartstat.x  tests_intpartstatntmp.par >/dev/null || exit 1

sed -e 's/ -3 / -4 /g' DATA_testn/testpart_proc000.control_part  >extract.dat 
sed -e 's/ -3 / -4 /g' DATA_testn_part/testpart_proc000.control_part  >extractpart.dat 

diff -q extract.dat extractpart.dat || { echo "code controle different particule !!!!" $group && exit 1 ; } 

for c in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 
do
    export PARTNO=PART$c

    grep "$PARTNO " DATA_testn/testpart_proc000.minmax_aei_part >extract.dat  || exit 1
    grep "$PARTNO " DATA_testn_part/testpart_proc000.minmax_aei_part >extractpart.dat || exit 1

    diff -w -q extract.dat  extractpart.dat || { echo "erreur min/max a,e,i n particule !!!!" $PARTNO && exit 1 ; }   
    rm extract.dat extractpart.dat

    grep "$PARTNO " DATA_testn/testpart_proc000.minmax_diffalp_part   >extract.dat  || exit 1
    grep "$PARTNO " DATA_testn_part/testpart_proc000.minmax_diffalp_part >extractpart.dat || exit 1

    diff -w -q extract.dat  extractpart.dat || { echo "erreur min/max alp   n particule !!!!" $PARTNO && exit 1 ; }   
    rm extract.dat extractpart.dat

    grep "$PARTNO " DATA_testn/testpart_proc000.naf_diffalp_part   >extract.dat  || exit 1
    grep "$PARTNO " DATA_testn_part/testpart_proc000.naf_diffalp_part >extractpart.dat || exit 1

    diff -w -q extract.dat  extractpart.dat || wc -l extract.dat  extractpart.dat || { echo "erreur naf dif alp   n particule !!!!" $PARTNO && exit 1 ; }   
    rm extract.dat extractpart.dat

    grep "$PARTNO " DATA_testn/testpart_proc000.naf_alkhqp_part    >extract.dat 
    grep "$PARTNO " DATA_testn_part/testpart_proc000.naf_alkhqp_part >extractpart.dat 
    diff -w -q extract.dat  extractpart.dat  || { echo "erreur naf alkqp   n particule !!!!" $PARTNO && exit 1 ; } 

rm -f extractok.dat
(trip <<__EOF
try {
for j=1 to 45
{
 vnumR xpart;
 fmt="%s %g"$
 read("extractpart.dat", fmt, unused, (xpart, j+2) );
 vnumR xpla;
 fmt="%s %g"$
 read("extract.dat", fmt, unused, (xpla, j+2) );
 if (size(xpart)!=0) then { if (max(abs(xpla-xpart))!=0) then { j; error("erreur dans la colonne "); }; };
}; 
t=0,10$
write("extractok.dat", t);
} catch { msg "une erreur a eu lieu!!!";};
__EOF
) | grep -v BIENVENUE | grep -v Taper | grep -v apprecie | sed '/^$/d'
    ls extractok.dat >/dev/null ||  exit 1 
    rm -f extractok.dat

    echo "particule "$c".... ok"

done
done

echo
echo "controle distance planete n particule.... ok"
echo 
rm -r -f DATA_testn_part DATA_testn extract.dat extractpart.dat extractok.dat
exit 0


