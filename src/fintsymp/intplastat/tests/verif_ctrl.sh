#!/bin/bash
rm -r -f DATA_test 
mkdir DATA_test
echo -n "ctrl ener............." 
../intplastat.x  tests_intplastatctrlener1.par >/dev/null || exit 1
grep ' -5' DATA_test/ctrlener__proc000.control >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
echo "ok" 
rm -r -f DATA_test
mkdir DATA_test
echo -n "ctrl dist star min............." 
../intplastat.x  tests_intplastatctrldist1.par >/dev/null || exit 1
grep ' -6' DATA_test/ctrldistmin__proc000.control >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
echo "ok" 
rm -r -f DATA_test
mkdir DATA_test
echo -n "ctrl dist star max compute=1............." 
../intplastat.x  tests_intplastatctrldist2.par >/dev/null || exit 1
grep ' -7' DATA_test/ctrldistmax__proc000.control >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
echo "ok" 
rm -r -f DATA_test
mkdir DATA_test
echo -n "ctrl dist star helio max compute=2............." 
../intplastat.x  tests_intplastatctrldist3.par >/dev/null || exit 1
grep ' -7' DATA_test/ctrldistmax__proc000.control >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
grep 'P009295' DATA_test/ctrldistmax__proc000.ctrlstar_car >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
grep 'P009294' DATA_test/ctrldistmax__proc000.ctrlstar_car >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
echo "ok" 
rm -r -f DATA_test
mkdir DATA_test
echo -n "ctrl dist star helio max compute=3............." 
../intplastat.x  tests_intplastatctrldist4.par >/dev/null || exit 1
grep ' -7' DATA_test/ctrldistmax__proc000.control >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
grep 'P009295' DATA_test/ctrldistmax__proc000.ctrlstar_ell >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
grep 'P009294' DATA_test/ctrldistmax__proc000.ctrlstar_ell >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
echo "ok" 
rm -r -f DATA_test
mkdir DATA_test
echo -n "ctrl dist star jacobi max compute=2............." 
../intplastat.x  tests_intplastatctrldist5.par >/dev/null || exit 1
grep ' -7' DATA_test/ctrldistmax__proc000.control >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
grep 'P009295' DATA_test/ctrldistmax__proc000.ctrlstar_car >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
grep 'P009294' DATA_test/ctrldistmax__proc000.ctrlstar_car >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
echo "ok" 
rm -r -f DATA_test
mkdir DATA_test
echo -n "ctrl dist star jacobi max compute=3............." 
../intplastat.x  tests_intplastatctrldist6.par >/dev/null || exit 1
grep ' -7' DATA_test/ctrldistmax__proc000.control >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
grep 'P009295' DATA_test/ctrldistmax__proc000.ctrlstar_ell >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
grep 'P009294' DATA_test/ctrldistmax__proc000.ctrlstar_ell >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
echo "ok" 

# cas de la distance minimale entre planetes
rm -r -f DATA_test
mkdir DATA_test
echo -n "ctrl dist pla min compute=1............." 
../intplastat.x  tests_intplastatctrldistpla1.par >/dev/null || exit 1
grep ' -9' DATA_test/ctrldistmin__proc000.control >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
grep '115' DATA_test/ctrldistmin__proc000.control >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
echo "ok" 
rm -r -f DATA_test
mkdir DATA_test
echo -n "ctrl dist pla helio compute=2............." 
../intplastat.x  tests_intplastatctrldistpla2.par >/dev/null || exit 1
grep ' -9' DATA_test/ctrldistmin__proc000.control >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
grep 'P009294' DATA_test/ctrldistmin__proc000.ctrlpla_car >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
echo "ok" 
rm -r -f DATA_test
mkdir DATA_test
echo -n "ctrl dist pla helio compute=3............." 
../intplastat.x  tests_intplastatctrldistpla3.par >/dev/null || exit 1
grep ' -9' DATA_test/ctrldistmin__proc000.control >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
grep 'P009294' DATA_test/ctrldistmin__proc000.ctrlpla_ell >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
echo "ok" 
rm -r -f DATA_test
mkdir DATA_test
echo -n "ctrl dist pla jacobi compute=2............." 
../intplastat.x  tests_intplastatctrldistpla4.par >/dev/null || exit 1
grep ' -9' DATA_test/ctrldistmin__proc000.control >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
grep 'P009294' DATA_test/ctrldistmin__proc000.ctrlpla_car >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
echo "ok" 
rm -r -f DATA_test
mkdir DATA_test
echo -n "ctrl dist pla jacobi compute=3............." 
../intplastat.x  tests_intplastatctrldistpla5.par >/dev/null || exit 1
grep ' -9' DATA_test/ctrldistmin__proc000.control >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
grep 'P009294' DATA_test/ctrldistmin__proc000.ctrlpla_ell >/dev/null || ( echo "FAIL" && exit 77 ); [ "$?" -eq 77 ] && exit 1
echo "ok" 
rm -r -f DATA_test

exit 0


