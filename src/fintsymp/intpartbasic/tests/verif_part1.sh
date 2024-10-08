#!/bin/bash
for outvell in 1 2 3 4 6 7 8
do
echo 'out_ell=' $outvell
sed -e "s/out_ell  = 4/out_ell = $outvell/" intparttest1.par >intparttest2.par
sed -e "s/out_ell  = 4/out_ell = $outvell/" intplabasictest1.par >intplabasictest2.par
rm -r -f DATA_testref DATA_test1
mkdir DATA_testref
../../intplabasic/intplabasic.x  intplabasictest2.par >/dev/null || exit 1
mkdir DATA_test1
../intpartbasic.x intparttest2.par >/dev/null || exit 1
diff DATA_test1/testpart_N0001_ABAH4.int DATA_testref/testpla_N0001_ABAH4.int || { echo "energie.....................FAILED!" && exit 1 ; }
echo "energie.....................ok"

echo 0 >testok
(trip<<__EOF
myerror=0$
vnumR t1, xyzpla[1:24];
read("DATA_testref/testpla_N0001_ABAH4.car",t1, xyzpla);
vnumR xyz1[1:6];
xyz1[1:6]=xyzpla[19:24];
fmt="%s "+7*"%g "$
vnumR t2, xyz2[1:6];
read("DATA_test1/testpart_PART001_ABAH4_part.car", fmt, ts, t2,xyz2);

for j=1 to 6 {
 if (max(abs(xyz1[j]-xyz2[j]))>2E-16) then {
    j;
    myerror=1$
    !"echo 1>testok";
    max(abs(xyz1[j]-xyz2[j]));
    msg "-----------------------------------------";
    msg "-----------------------------------------";
    msg "erreur dans les coordonnees cartesiennes de la particule intparttest1.par";
    msg "-----------------------------------------";
    error("-----------------------------------------");
 };

};
if (myerror==0) then 
{
    msg "coordonnes cartesiennes pour 1 particule .....................ok";
};
vnumR t1, xyzpla[1:24];
read("DATA_testref/testpla_N0001_ABAH4.ell",t1, xyzpla);
vnumR xyz1[1:6];
xyz1[1:6]=xyzpla[19:24];
fmt="%s "+7*"%g "$
vnumR t2, xyz2[1:6];
read("DATA_test1/testpart_PART001_ABAH4_part.ell", fmt, ts, t2,xyz2);

for j=1 to 6 {
 if (max(abs(xyz1[j]-xyz2[j]))>2E-16) then {
    myerror=1$
    !"echo 1>testok";
    j;
    msg "-----------------------------------------";
    msg "-----------------------------------------";
    msg "erreur dans les coordonnees elliptiques de la particule intparttest1.par";
    max(abs(xyz1[j]-xyz2[j]));
    msg "-----------------------------------------";
    error("-----------------------------------------");
 };

};
if (myerror==0) then 
{
    msg "coordonnes elliptiques pour 1 particule .....................ok";
};
__EOF
) | grep -v BIENVENUE | grep -v Taper | grep -v apprecie 

grep 0 testok >/dev/null || { echo "erreur dans les coordonnees ... de la particule intparttest1.par !!!" && exit 1 ; }

rm -r -f DATA_testref DATA_test1 testok
done

exit 0


