#!/bin/bash
for outvell in 1 2 3 4 6 7 8
do
echo 'out_ell=' $outvell
sed -e "s/out_ell  = 4/out_ell = $outvell/" intparttestn.par >intparttestn2.par
sed -e "s/out_ell  = 4/out_ell = $outvell/" intplabasictestn.par >intplabasictestn2.par
rm -r -f DATA_testnref DATA_testn
mkdir DATA_testnref
../../intplabasic/intplabasic.x  intplabasictestn2.par >/dev/null || exit 1
mkdir DATA_testn
../intpartbasic.x intparttestn2.par >/dev/null || exit 1
diff DATA_testn/testpart_N0001_ABAH4.int DATA_testnref/testpla_N1_ABAH4.int || { echo "energie.....................FAILED!" && exit 1 ; }
echo "energie.....................ok"

echo 0 >testok

(trip<<__EOF
msg "verification par block de 1.................... demarrage";
myerror=0$
for vipart=0 to 10 {
ipart=str(vipart)$
vnumR t1, xyzpla[1:24];
read("DATA_testnref/testpla_N"+ipart+"_ABAH4.car",t1, xyzpla);
vnumR xyz1[1:6];
xyz1[1:6]=xyzpla[19:24];
fmt="%s "+7*"%g "$
vnumR t2, xyz2[1:6];
read("DATA_testn/testpart_"+ipart+"_ABAH4_part.car", fmt, ts, t2,xyz2);

for j=1 to 6 {
 if (max(abs(xyz1[j]-xyz2[j]))>2E-16) then {
    j;
    myerror=1$
    !"echo 1>testok";
    msg "-----------------------------------------";
    msg "-----------------------------------------";
    msg "erreur dans les coordonnees cartesiennes de la particule  "+ipart;
    max(abs(xyz1[j]-xyz2[j]));
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
read("DATA_testn/testpart_"+ipart+"_ABAH4_part.ell", fmt, ts, t2,xyz2);

for j=1 to 6 {
 if (max(abs(xyz1[j]-xyz2[j]))>2e-16) then {
    j;
    myerror=1$
    !"echo 1>testok";
    msg "-----------------------------------------";
    msg "-----------------------------------------";
    msg "erreur dans les coordonnees elliptiques de la particule "+ipart;
    max(abs(xyz1[j]-xyz2[j]));
    msg "-----------------------------------------";
    error("-----------------------------------------");
 };

};
};
if (myerror==0) then 
{
msg "coordonnes cartesiennes pour particule par block de 1....................ok";
msg "coordonnes elliptiques pour particule par block de 1....................ok";
};
__EOF
) | grep -v BIENVENUE | grep -v Taper | grep -v apprecie 

grep 0 testok >/dev/null || { echo "erreur dans les coordonnees elliptiques de la particule !!!" && exit 1 ; }

rm -r -f DATA_testnref DATA_testn
done


for outvell in 1 2 3 4 6 7 8
do
echo 'out_ell=' $outvell
for partblock in 1 2 3 4 5 6  7 8 9 10
do
echo $partblock 'particules par bloc - out_ell=', $outvell
sed -e "s/out_ell  = 4/out_ell = $outvell/" intparttestn.par | sed -e "s/part_blocksize = 1/part_blocksize = $partblock/" >intparttestn2.par
sed -e "s/out_ell  = 4/out_ell = $outvell/" intplabasictestn.par >intplabasictestn2.par
rm -r -f DATA_testnref DATA_testn
mkdir DATA_testnref
../../intplabasic/intplabasic.x  intplabasictestn2.par >/dev/null || exit 1
mkdir DATA_testn
../intpartbasic.x intparttestn2.par >/dev/null || exit 1
diff DATA_testn/testpart_N0001_ABAH4.int DATA_testnref/testpla_N1_ABAH4.int || { echo "energie.....................FAILED!" && exit 1 ; }
echo "energie.....................ok"
cat DATA_testn/testpart_*_ABAH4_part.car >DATA_testn/testpart_PART_full_part.car
cat DATA_testn/testpart_*_ABAH4_part.ell >DATA_testn/testpart_PART_full_part.ell

(trip<<__EOF
blocksize=${partblock}$
myerror=0$
msg ("verification par block de %d.................... demarrage\n", blocksize);
vnumR tsall, t2all, xyz2all[1:6];
read("DATA_testn/testpart_PART_full_part.car",tsall, t2all,xyz2all);
vnumR tsell, t2ell, xyz2ell[1:6];
read("DATA_testn/testpart_PART_full_part.ell", tsell, t2ell,xyz2ell);
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
 if (max(abs(xyz1[j]-xyz2[j]))>2E-16) then {
    j;
    myerror=1$
    !"echo 1>testok";
    msg "-----------------------------------------";
    msg "-----------------------------------------";
    msg "erreur dans les coordonnees cartesiennes de la particule  "+ipart;
    max(abs(xyz1[j]-xyz2[j]));
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
 if (max(abs(xyz1[j]-xyz2[j]))>2E-16) then {
    j;
    !"echo 1>testok";
    myerror=1$
    msg "-----------------------------------------";
    msg "-----------------------------------------";
    msg "erreur dans les coordonnees elliptiques de la particule "+ipart;
    max(abs(xyz1[j]-xyz2[j]));
    msg "-----------------------------------------";
    error("-----------------------------------------");
 };

};
};
if (myerror==0) then 
{
msg ("coordonnes cartesiennes pour particule par block de %d.................. ok\n", blocksize);
msg ("coordonnes elliptiques pour particule par block de %d.................. ok\n", blocksize);
};
__EOF
) | grep -v BIENVENUE | grep -v Taper | grep -v apprecie 

grep 0 testok >/dev/null || { echo "erreur dans les coordonnees elliptiques de la particule !!!" && exit 1 ; }

rm -r -f DATA_testnref DATA_testn 
done
done
rm  -f  testok

exit 0


