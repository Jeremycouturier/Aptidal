#!/bin/bash
rm -r -f DATA_test2 DATA_test testok
mkdir DATA_test
../../intplastat/intplastat.x  tests_intplastatmaxaei1h.par >/dev/null || exit 1
mkdir DATA_test2
cp tests_filepartci.dat tests_filepartci1.dat
../intpartstat.x  tests_intpartstatnaf1h.par >/dev/null || exit 1

sed -i -e 's/SYS001/PART1/g'  DATA_test/test2_proc000.naf_diffalp
echo 0 >testok

(trip <<__EOF
n=16$ vnumR w1[1:n], w2[1:n];
fmt="%s "+n*"%g"$
read("DATA_test/test2_proc000.naf_diffalp", fmt, s, w1);
read("DATA_test2/test_proc000.naf_diffalp_part", fmt, s, w2);
for j=1 to n {
 if (max(abs(w1[j]-w2[j]))>1) then {  max(abs(w1[j]-w2[j])/abs(w2[j]));
 msg("probleme dans colonne %d\n", j); !"echo 1 >testok"; };
};	
__EOF
) | grep -v BIENVENUE | grep -v Taper | grep -v apprecie | sed '/^$/d'
grep 0 testok >/dev/null || { echo "eerreur naf diff angle 1 particule!!!!" && exit 1 ; }

#diff -w -q DATA_test/test2_proc000.naf_diffalp DATA_test2/test_proc000.naf_diffalp_part || { echo "erreur naf diff angle 1 particule !!!!" && exit 1 ; }

echo "naf diff angle  1 particule.... ok"

rm -r -f DATA_test2 DATA_test testok
exit 0


