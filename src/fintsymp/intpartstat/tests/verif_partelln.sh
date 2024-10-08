#!/bin/bash
for outvell in 1 2 3 4 6 7 8
do
echo 'out_ell=' $outvell
for partblock in 1 2 3 4 5 6  7 8 9 10
do
echo $partblock 'particules par bloc - out_ell=' $outvell
sed -e "s/out_ell_part  = 4/out_ell_part = $outvell/" tests_intpartstatell.par | sed -e "s/part_blocksize = 1/part_blocksize = $partblock/" >intparttestn2.par
sed -e "s/out_ell  = 4/out_ell = $outvell/" intplabasictestelln2.par >intplabasictestn2.par
rm -r -f DATA_testnref DATA_test2
mkdir DATA_testnref
../../intplabasic/intplabasic.x  intplabasictestn2.par >/dev/null || exit 1
mkdir DATA_test2
../intpartstat.x intparttestn2.par >/dev/null || exit 1
diff DATA_test2/test_proc000.int DATA_testnref/testpla_N0_ABAH4.int || { echo "energie.....................FAILED!" && exit 1 ; }
echo "energie.....................ok"


rm -f extractok.dat
(trip <<__EOF
try { 
blocksize=${partblock}$
msg ("verification par block de %d.................... demarrage\n", blocksize);
vnumR tsall, t2all, xyz2all[1:6];
read("DATA_test2/test_proc000.car_part",tsall, t2all,xyz2all);
vnumR tsell, t2ell, xyz2ell[1:6];
read("DATA_test2/test_proc000.ell_part", tsell, t2ell,xyz2ell);
for vipart=0 to 10 {
ipart=str(vipart)$
vnumR t1, xyzpla[1:24];
read("DATA_testnref/testpla_N"+ipart+"_ABAH4.car",t1, xyzpla);
vnumR xyz1[1:6];
xyz1[1:6]=xyzpla[19:24];
vnumR t2, xyz2[1:6];
t2=select(tsall==vipart, t2all)$
for j=1 to 6 { xyz2[j]=select(tsall==vipart, xyz2all[j])$ };

for j=1 to 6 {
 if (max(abs(xyz1[j]-xyz2[j]))!=0) then {
    j;
    msg "-----------------------------------------";
    msg "-----------------------------------------";
    msg "erreur dans les coordonnees cartesiennes de la particule  "+ipart;
    msg "-----------------------------------------";
    error("-----------------------------------------");
 };

};

vnumR t1, xyzpla[1:24];
read("DATA_testnref/testpla_N"+ipart+"_ABAH4.ell",t1, xyzpla);
vnumR xyz1[1:6];
xyz1[1:6]=xyzpla[19:24];
fmt="%s "+7*"%g "$
vnumR t2, xyz2[1:6];
t2=select(tsell==vipart, t2ell)$
for j=1 to 6 { xyz2[j]=select(tsell==vipart, xyz2ell[j])$ };

for j=1 to 6 {
 if (max(abs(xyz1[j]-xyz2[j]))!=0) then {
    j;
    msg "-----------------------------------------";
    msg "-----------------------------------------";
    msg "erreur dans les coordonnees elliptiques de la particule "+ipart;
    max(abs(xyz1[j]-xyz2[j]));
    msg "-----------------------------------------";
    error("-----------------------------------------");
 };

};
};
msg ("coordonnes cartesiennes pour particule par block de %d.................. ok\n", blocksize);
msg ("coordonnes elliptiques pour particule par block de %d.................. ok\n", blocksize);
t=0,10$
write("extract.dat", t);
} catch { msg "une erreur a eu lieu!!!";};
__EOF
) | grep -v BIENVENUE | grep -v Taper | grep -v apprecie | sed '/^$/d'
ls extract.dat >/dev/null || { echo "erreur particule par block !!!!" && exit 1 ; }
rm -f extract.dat

rm -r -f DATA_testnref DATA_test2 intplabasictestn2.par intparttestn2.par
done
done

exit 0


