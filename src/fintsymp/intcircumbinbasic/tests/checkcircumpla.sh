#! /bin/bash
rm -r -f ./DATA verif.dat
mkdir -p ./DATA
../intcircumbinbasic.x checkcircum_file0.par >/dev/null
../intcircumbinbasic.x checkcircum_file2.par >/dev/null
../../intplabasic/intplabasic.x checkbasic_file0.par >/dev/null
../../intplabasic/intplabasic.x checkbasic_file2.par >/dev/null

../intcircumbinbasic.x checkcircum_fileinv0.par >/dev/null
../intcircumbinbasic.x checkcircum_fileinv2.par >/dev/null
../../intplabasic/intplabasic.x checkbasic_fileinv0.par >/dev/null
../../intplabasic/intplabasic.x checkbasic_fileinv2.par >/dev/null

echo "verification vs planetaire heliocentrique : 1 planete dans le repere original....."
echo 0 >verif.dat
trip >/dev/null <<__EOF
vnumR win[1:12], wout[1:12], t;
read("DATA/verifplabasic0_1_ABAH4.car", t, win );
read("DATA/verifplacircum0_1_ABAH4.car",t, wout );
t=1,1$
if (max(abs(win-wout))<=1E-11) then { write(verif.dat, t); };
__EOF
grep 1 verif.dat >/dev/null ||  { echo "erreur !!!" && exit 1 ; }
echo "ok"

echo "verification vs planetaire heliocentrique : 2 planetes dans le repere original....."
echo 0 >verif.dat
trip >/dev/null <<__EOF
vnumR win[1:18], wout[1:18], t;
read("DATA/verifplabasic2_1_ABAH4.car",t, win );
read("DATA/verifplacircum2_1_ABAH4.car",t, wout );
t=1,1$
if (max(abs(win-wout))<=1E-11) then { write(verif.dat, t); };
__EOF
grep 1 verif.dat >/dev/null ||  { echo "erreur !!!" && exit 1 ; }
echo "ok"

echo "verification vs planetaire heliocentrique des integrales premieres : 2 planetes dans le repere original....."
echo 0 >verif.dat
trip >/dev/null <<__EOF
vnumR win[1:5], wout[1:5], t, m[1:2];
read("DATA/verifplabasic2_1_ABAH4.int",win );
read("DATA/verifplacircum2_1_ABAH4.int", wout );
read("checkbasic_file2.dat",[1:1],(m,3) );
t=1,1$
ms=m[1][1];
for j=2 to 5 {
if (max(abs(win[j]-wout[j]/ms))>=1E-14) then { t=t*0; };
};
write(verif.dat, t);
__EOF
grep 1 verif.dat >/dev/null ||  { echo "erreur !!!" && exit 1 ; }
echo "ok"

echo "verification vs planetaire heliocentrique : 1 planete dans le repere invariant....."
echo 0 >verif.dat
trip >/dev/null <<__EOF
vnumR win[1:12], wout[1:12], t;
read("DATA/verifplabasicinv0_1_ABAH4.car", t, win );
read("DATA/verifplacircuminv0_1_ABAH4.car",t, wout );
t=1,1$
if (max(abs(win-wout))<=1E-11) then { write(verif.dat, t); };
__EOF
grep 1 verif.dat >/dev/null ||  { echo "erreur !!!" && exit 1 ; }
echo "ok"

echo "verification vs planetaire heliocentrique : 2 planetes  dans le repere invariant....."
echo 0 >verif.dat
trip >/dev/null <<__EOF
vnumR win[1:18], wout[1:18], t;
read("DATA/verifplabasicinv2_1_ABAH4.car",t, win );
read("DATA/verifplacircuminv2_1_ABAH4.car",t, wout );
t=1,1$
if (max(abs(win-wout))<=1E-11) then { write(verif.dat, t); };
__EOF
grep 1 verif.dat >/dev/null ||  { echo "erreur !!!" && exit 1 ; }
echo "ok"
echo "fin..."
exit 0
