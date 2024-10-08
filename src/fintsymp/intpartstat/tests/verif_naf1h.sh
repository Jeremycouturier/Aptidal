#!/bin/bash
rm -r -f DATA_test2 DATA_test
cp tests_fileci.dat tests_fileci1.dat
mkdir DATA_test
../../intplastat/intplastat.x  tests_intplastatnaf1h.par >/dev/null || exit 1
mkdir DATA_test2
cp tests_filepartci.dat tests_filepartci1.dat
../intpartstat.x  tests_intpartstatnaf1h.par >/dev/null || exit 1

diff DATA_test2/test_proc000.naf_alkhqp DATA_test/test_proc000.naf_alkhqp || { echo "erreur analyse en frequence des planetes !!!!" && exit 1 ; }

echo "analyse en frequence des planetes.... ok"

rm -r -f DATA_test2 DATA_test
exit 0


