#!/bin/bash
rm -r -f DATA_test2 DATA_test
cp tests_fileci.dat tests_fileci1.dat
mkdir DATA_test
../intplastat.x  tests_intplastatminmaxae2j.par >/dev/null || exit 1

#verification des min/max (a_p(1)-a_p(2))**2, ...
rm -f extractok.dat
(trip<<__EOF
try {
npla=2$
fminmax_ae2="./DATA_test/testj_proc000.minmax_diffae2"$
fmt="%s %g "+9*" %g"$
vnumR t_minmax_ae2,minmax_ae2[1:9]$
read(fminmax_ae2, fmt, rad_minmax_ae2, t_minmax_ae2, minmax_ae2);
for ci=1 to 3 {
rad=str("N000%d", ci)$
fell="./DATA_test/testj_"+rad+".ell"$
vnumR t, ell[1:npla*6];
read(fell, t, ell);

id=1,size(t_minmax_ae2)$
q=dimtovnumR(select(rad_minmax_ae2==rad,vnumtodim(id)))$

vnumR vt[1:5]([1:size(q)])$
vnumR vmin[1:5]([1:size(q)])$
vnumR vmax[1:5]([1:size(q)])$
vnumR vmoy[1:5]([1:size(q)])$
t0=0$
for k=1 to size(q)
{
 t1 = t_minmax_ae2[q[k]]$
 for j=1 to 1 {
  v1 = select((t0<=t) && (t<=t1), ell[6*0+j])$
  v2 = select((t0<=t) && (t<=t1), ell[6*1+j])$
  v=(v1-v2)**2$
  vt[j][k] = t1$
  vmin[j][k]=min(v)$
  vmax[j][k]=max(v)$
  vmoy[j][k]=sum(v)/size(v)$
 };
 t0=t1$
};

 for j=1 to 1 {
    diffmin=max(abs(minmax_ae2[3*(j-1)+1][q]-vmin[j]))$
    diffmax=max(abs(minmax_ae2[3*(j-1)+3][q]-vmax[j]))$
    diffmoy=max(abs(minmax_ae2[3*(j-1)+2][q]-vmoy[j]))$
    if ((diffmin>1E-12) || (diffmax>1E-12) || (diffmoy>1E-12)) then 
    {
        error("erreur min/max de (a_p(1)-a_p(2))**2 - jacobi");
    };
 };

t0=0$
for k=1 to size(q)
{
  t1 = t_minmax_ae2[q[k]]$
  v1 = select((t0<=t) && (t<=t1), sqrt(ell[6*0+3]**2+ell[6*0+4]**2))$
  v2 = select((t0<=t) && (t<=t1), sqrt(ell[6*1+3]**2+ell[6*1+4]**2))$
  v=(v1-v2)**2$
  j=2$
  vt[j][k] = t1$
  vmin[j][k]=min(v)$
  vmax[j][k]=max(v)$
  vmoy[j][k]=sum(v)/size(v)$
 t0=t1$
};

 for j=2 to 2 {
    diffmin=max(abs(minmax_ae2[3*(j-1)+1][q]-vmin[j]))$
    diffmax=max(abs(minmax_ae2[3*(j-1)+3][q]-vmax[j]))$
    diffmoy=max(abs(minmax_ae2[3*(j-1)+2][q]-vmoy[j]))$
    if ((diffmin>1E-12) || (diffmax>1E-12) || (diffmoy>1E-12)) then 
    {
        error("erreur min/max de (e_p(1)-e_p(2))**2 - jacobi");
    };
 };


t0=0$
for k=1 to size(q)
{
  t1 = t_minmax_ae2[q[k]]$
  v3 = select((t0<=t) && (t<=t1), ell[6*0+1])$
  v4 = select((t0<=t) && (t<=t1), ell[6*1+1])$
  v1 = select((t0<=t) && (t<=t1), sqrt(ell[6*0+3]**2+ell[6*0+4]**2))$
  v2 = select((t0<=t) && (t<=t1), sqrt(ell[6*1+3]**2+ell[6*1+4]**2))$
  v=(v3-v4)**2+(v1-v2)**2$
  j=3$
  vt[j][k] = t1$
  vmin[j][k]=min(v)$
  vmax[j][k]=max(v)$
  vmoy[j][k]=sum(v)/size(v)$
 t0=t1$
};

 for j=3 to 3 {
    diffmin=max(abs(minmax_ae2[3*(j-1)+1][q]-vmin[j]))$
    diffmax=max(abs(minmax_ae2[3*(j-1)+3][q]-vmax[j]))$
    diffmoy=max(abs(minmax_ae2[3*(j-1)+2][q]-vmoy[j]))$
    if ((diffmin>1E-12) || (diffmax>1E-12) || (diffmoy>1E-12)) then 
    {
        error("erreur min/max de (a_p(1)-a_p(2))**2+(e_p(1)-e_p(2))**2 - jacobi");
    };
 };

};
msg "min/max/moy (jacobi) (a_p(1)-a_p(2))**2, ... .....................ok";
t=0,10$
write("extractok.dat", t);
}
catch {
msg "-----------------------------------------";
msg "-----------------------------------------";
msg "erreur dans la determination du min/max de (a_p(1)-a_p(2))**2";
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


