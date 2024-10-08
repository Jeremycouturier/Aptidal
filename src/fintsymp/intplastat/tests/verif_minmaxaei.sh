#!/bin/bash
rm -r -f DATA_test2 DATA_test
cp tests_fileci.dat tests_fileci1.dat
mkdir DATA_test
../intplastat.x  tests_intplastath.par >/dev/null || exit 1
mv DATA_test DATA_test2
mkdir DATA_test
cp DATA_test2/test_proc000.ci tests_fileci1.dat
../intplastat.x  tests_intplastath.par >/dev/null || exit 1
rm -f DATA_test2/test.par DATA_test/test.par
diff -r -q DATA_test DATA_test2 || exit 1

#verification des min/max a,e,i
rm -f extractok.dat
(trip<<__EOF
npla=4$
try {
fminmax_aei="./DATA_test/test_proc000.minmax_aei"$
fmt="%s %g "+9*npla*" %g"$
vnumR t_minmax_aei,minmax_aei[1:9*npla]$
read(fminmax_aei, fmt, rad_minmax_aei, t_minmax_aei, minmax_aei);
for ci=1 to 3 {
rad=str("N000%d", ci)$
fell="./DATA_test/test_"+rad+".ell"$
vnumR t, ell[1:npla*6];
read(fell, t, ell);

id=1,size(t_minmax_aei)$
q=dimtovnumR(select(rad_minmax_aei==rad,vnumtodim(id)))$

vnumR vt[1:npla,1:3]([1:size(q)])$
vnumR vmin[1:npla,1:3]([1:size(q)])$
vnumR vmax[1:npla,1:3]([1:size(q)])$
vnumR vmoy[1:npla,1:3]([1:size(q)])$
for kpla=0 to npla-1
{
t0=0$
kpla_cur = kpla$
for k=1 to size(q)
{
 t1 = t_minmax_aei[q[k]]$
 for j=1 to 3 {
  v = select((t0<=t) && (t<=t1), ell[6*kpla+j])$
  vt[kpla+1,j][k] = t1$
  vmin[kpla+1,j][k]=min(v)$
  vmax[kpla+1,j][k]=max(v)$
  vmoy[kpla+1,j][k]=sum(v)/size(v)$
 };
 t0=t1$
};

 for j=1 to 3 {
    diffmin=max(abs(minmax_aei[9*kpla+3*(j-1)+1][q]-vmin[kpla+1,j]))$
    diffmax=max(abs(minmax_aei[9*kpla+3*(j-1)+3][q]-vmax[kpla+1,j]))$
    diffmoy=max(abs(minmax_aei[9*kpla+3*(j-1)+2][q]-vmoy[kpla+1,j]))$
    if ((diffmin>1E-12) || (diffmax>1E-12) || (diffmoy>1E-12)) then 
    {
        error("erreur min/max de a,e,I");
    };
 };
};
};
msg "min/max/moy a,e,I.....................ok";
t=0,10$
write("extractok.dat", t);
}
catch {
msg "-----------------------------------------";
msg "-----------------------------------------";
msg "erreur dans la determination du min/max de a,e,I";
kpla_cur;
diffmin;
diffmax;
diffmoy;
msg "-----------------------------------------";
error("-----------------------------------------");
};
__EOF
) | grep -v BIENVENUE | grep -v Taper | grep -v apprecie | sed '/^$/d'
    ls extractok.dat >/dev/null ||  exit 1 
    rm -f extractok.dat
rm -r -f DATA_test2 DATA_test
exit 0


