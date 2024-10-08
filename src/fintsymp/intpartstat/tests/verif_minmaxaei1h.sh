#!/bin/bash
rm -r -f DATA_test2 DATA_test testok
mkdir DATA_test
../../intplastat/intplastat.x  tests_intplastatmaxaei1h.par >/dev/null || exit 1
mkdir DATA_test2
cp tests_filepartci.dat tests_filepartci1.dat
../intpartstat.x  tests_intpartstatnaf1h.par >/dev/null || exit 1

awk '{ print "PART1 "$2" "$30" "$31" "$32" "$33" "$34" "$35" "$36" "$37" "$38 }' DATA_test/test2_proc000.minmax_aei >extract.dat

echo 0 >testok
trip  >/dev/null  <<__EOF
vnumR w1[1:10], w2[1:10];
fmt="%s "+10*"%g";
read("extract.dat", fmt, s, w1);
read("DATA_test2/test_proc000.minmax_aei_part", fmt, s, w2);
for j=1 to 10 {
 if (max(abs(w1[j]-w2[j])/abs(w2[j]))>2E-14) then {  max(abs(w1[j]-w2[j])/abs(w2[j]));
 msg("probleme dans colonne %d\n", j); !"echo 1 >testok"; };
};	
__EOF
grep 0 testok >/dev/null  || { echo "erreur min/max a,e,i 1 particule !!!!" && exit 1 ; }

# diff -w -q extract.dat DATA_test2/test_proc000.minmax_aei_part  || { echo "erreur min/max a,e,i 1 particule !!!!" && exit 1 ; }

echo "min/max a,e,i 1 particule.... ok"

rm -r -f DATA_test2 DATA_test extract.dat
exit 0


