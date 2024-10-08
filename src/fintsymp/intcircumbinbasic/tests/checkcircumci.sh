#! /bin/bash
rm -r -f ./DATA verif.dat
mkdir -p ./DATA
../intcircumbinbasic.x checkcircumci_file13.par >/dev/null
../intcircumbinbasic.x checkcircumci_file11.par >/dev/null
../intcircumbinbasic.x checkcircumci_file12.par >/dev/null
../intcircumbinbasic.x checkcircumci_file14.par >/dev/null
../intcircumbinbasic.x checkcircumout_file14.par >/dev/null
../intcircumbinbasic.x checkcircumout_file13.par >/dev/null

echo "element elliptique : type 11....."
echo 0 >verif.dat
trip >/dev/null <<__EOF
vnumR win[1:12], wout[1:12], t;
read("DATA/verifci11__1_ABAH4.ell",[1:1], t, wout );
read("checkcircumci_file11.dat",[1:1], (win, 7) );
t=1,1$
if (max(abs(win-wout))<=1E-14) then { write(verif.dat, t); };
__EOF
grep 1 verif.dat >/dev/null ||  { echo "erreur !!!" && exit 1 ; }
echo "ok"
echo "element elliptique : type 13....."
echo 0 >verif.dat
trip >/dev/null <<__EOF
vnumR win[1:12], wout[1:12], t;
read("DATA/verifci13__1_ABAH4.ell",[1:1], t, wout );
read("checkcircumci_file13.dat",[1:1], (win, 7) );
t=1,1;
if (max(abs(win-wout))<=1E-14) then { write(verif.dat, t); };
__EOF
grep 1 verif.dat >/dev/null ||  { echo "erreur !!!" && exit 1 ; }
echo "ok"
echo "element elliptique : type 12....."
echo 0 >verif.dat
trip >/dev/null <<__EOF
vnumR win[1:12], wout[1:12], t;
read("DATA/verifci12__1_ABAH4.ell",[1:1], t, wout );
read("checkcircumci_file12.dat",[1:1], (win, 7) );
t=1,1$
if (max(abs(win-wout))<=1E-14) then { write(verif.dat, t); };
__EOF
grep 1 verif.dat >/dev/null ||  { echo "erreur !!!" && exit 1 ; }
echo "ok"
echo "element elliptique : type 14....."
echo 0 >verif.dat
trip >/dev/null <<__EOF
vnumR win[1:12], wout[1:12], t;
read("DATA/verifci14__1_ABAH4.ell",[1:1], t, wout );
read("checkcircumci_file14.dat",[1:1], (win, 8) );
t=1,1;
if (max(abs(win-wout))<=1E-14) then { write(verif.dat, t); };
__EOF
grep 1 verif.dat >/dev/null ||  { echo "erreur !!!" && exit 1 ; }
echo "ok"
echo "element elliptique : type 13 dans invariant....."
../intcircumbinbasic.x checkcircumci_fileinv13.par >/dev/null
echo 0 >verif.dat
trip >/dev/null <<__EOF
vnumR win[1:12], wout[1:12], t;
read("DATA/verifciinv13__1_ABAH4.ell",[1:1], t, wout );
read("checkcircumci_fileinv13.dat",[1:1], (win, 7) );
t=1,1;
if (max(abs(win-wout))<=1E-14) then { write(verif.dat, t); };
__EOF
grep 1 verif.dat >/dev/null ||  { echo "erreur !!!" && exit 1 ; }
echo "ok"
echo "element elliptique : type 14 ---> type 13....."
echo 0 >verif.dat
trip >/dev/null <<__EOF
vnumR win[1:12], wout[1:12], t;
read("DATA/verifci14to13__1_ABAH4.ell",[1:1], t, wout );
read("checkcircumout_file14.dat",[1:1], (win, 8) );
t=1,1;
if (max(abs(win-wout))<=1E-14) then { write(verif.dat, t); };
__EOF
echo "ok"
echo "element elliptique : type 13 ---> type 14....."
echo 0 >verif.dat
trip >/dev/null <<__EOF
vnumR win[1:12], wout[1:12], t;
read("DATA/verifci13to14__1_ABAH4.ell",[1:1], t, wout );
read("checkcircumci_file14.dat",[1:1], (win, 8) );
t=1,1;
if (max(abs(win-wout))<=1E-14) then { write(verif.dat, t); };
__EOF
echo "ok"
echo "fin..."
exit 0
