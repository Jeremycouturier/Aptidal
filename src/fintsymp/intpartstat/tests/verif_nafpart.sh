#!/bin/bash
rm -r -f DATA_test2 DATA_test
mkdir DATA_test
../../intplastat/intplastat.x  tests_intplastatmaxaei1h.par >/dev/null || exit 1
mkdir DATA_test2
cp tests_filepartci.dat tests_filepartci1.dat
../intpartstat.x  tests_intpartstatnaf1h.par >/dev/null || exit 1

rm -f extract.dat
(trip <<__EOF
try {
for j=1 to 45
{
 vnumR xpart;
 fmt="%s %g"$
 read("DATA_test2/test_proc000.naf_alkhqp_part", fmt, unused, (xpart, j+2) );
 vnumR xpla;
 fmt="%s %g"$
 read("DATA_test/test2_proc000.naf_alkhqp", fmt, unused, (xpla, j+2+45*3) );
 if (max(abs(xpla-xpart))>1E-10) then { j; max(abs(xpla-xpart)); error("erreur dans la colonne "); };
}; 
t=0,10$
write("extract.dat", t);
} catch { msg "une erreur a eu lieu!!!";};
__EOF
) | grep -v BIENVENUE | grep -v Taper | grep -v apprecie | sed '/^$/d'
ls extract.dat >/dev/null || { echo "erreur naf part 1 particule !!!!" && exit 1 ; }
rm -f extract.dat

echo "naf part  1 particule.... ok"
echo

rm -r -f DATA_test2 DATA_test extract.dat
exit 0


