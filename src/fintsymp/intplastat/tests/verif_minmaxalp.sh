#!/bin/bash
rm -r -f DATA_test2 DATA_test
cp tests_fileci.dat tests_fileci1.dat
mkdir DATA_test
../intplastat.x  tests_intplastatminmaxalp.par >/dev/null || exit 1

#verification des min/max a_p(1)-a_p(2), ... [-pi,pi] et [0,2*pi] 
rm -f extractok.dat
(trip<<__EOF
try {
npla=2$
fminmax_alp="./DATA_test/test_proc000.minmax_diffalp"$
fmt="%s %g "+15*" %g"$
vnumR t_minmax_alp,minmax_alp[1:15]$
read(fminmax_alp, fmt, rad_minmax_alp, t_minmax_alp, minmax_alp);
for ci=1 to 3 {
rad=str("N000%d", ci)$
fell="./DATA_test/test_"+rad+".ell"$
vnumR t, ell[1:npla*6];
read(fell, t, ell);

id=1,size(t_minmax_alp)$
q=dimtovnumR(select(rad_minmax_alp==rad,vnumtodim(id)))$

vnumR vt[1:5]([1:size(q)])$
vnumR vmin[1:5]([1:size(q)])$
vnumR vmax[1:5]([1:size(q)])$
vnumR vmoy[1:5]([1:size(q)])$
t0=0$
for k=1 to size(q)
{
 t1 = t_minmax_alp[q[k]]$
 for j=1 to 1 {
  v1 = select((t0<=t) && (t<=t1), ell[6*0+j])$
  v2 = select((t0<=t) && (t<=t1), ell[6*1+j])$
  v=v1-v2$
  vt[j][k] = t1$
  vmin[j][k]=min(v)$
  vmax[j][k]=max(v)$
  vmoy[j][k]=sum(v)/size(v)$
 };
 t0=t1$
};

 for j=1 to 1 {
    diffmin=max(abs(minmax_alp[3*(j-1)+1][q]-vmin[j]))$
    diffmax=max(abs(minmax_alp[3*(j-1)+3][q]-vmax[j]))$
    diffmoy=max(abs(minmax_alp[3*(j-1)+2][q]-vmoy[j]))$
    if ((diffmin>1E-12) || (diffmax>1E-12) || (diffmoy>1E-12)) then 
    {
        error("erreur min/max de a_p(1)-a_p(2) [-pi,pi] et [0,2*pi] ");
    };
 };

t0=0$
for k=1 to size(q)
{
  t1 = t_minmax_alp[q[k]]$
  v1 = select((t0<=t) && (t<=t1), ell[6*0+2])$
  v2 = select((t0<=t) && (t<=t1), ell[6*1+2])$
  v=mod(v1-v2, 2*pi)$
  v =?(v<0):v+2*pi:v$
  v_npi = mod(v1-v2, 2*pi)$
  v_npi =?(v_npi>pi):v_npi-2*pi:v_npi$
  v_npi =?(v_npi<-pi):v_npi+2*pi:v_npi$
  j=2$
  vt[j][k] = t1$
  vmin[j][k]=min(v)$
  vmax[j][k]=max(v)$
  vmoy[j][k]=sum(v)/size(v)$
  j=3$
  vmin[j][k]=min(v_npi)$
  vmax[j][k]=max(v_npi)$
  vmoy[j][k]=sum(v_npi)/size(v_npi)$
 t0=t1$
};

 for j=2 to 3 {
    diffmin=max(abs(minmax_alp[3*(j-1)+1][q]-vmin[j]))$
    diffmax=max(abs(minmax_alp[3*(j-1)+3][q]-vmax[j]))$
    diffmoy=max(abs(minmax_alp[3*(j-1)+2][q]-vmoy[j]))$
    if ((diffmin>1E-12) || (diffmax>1E-12) || (diffmoy>1E-12)) then 
    {
        error("erreur min/max de l_p(1)-l_p(2) [-pi,pi] et [0,2*pi] ");
    };
 };


t0=0$
for k=1 to size(q)
{
  t1 = t_minmax_alp[q[k]]$
  v1 = select((t0<=t) && (t<=t1), atan2(ell[6*0+4],ell[6*0+3]))$
  v2 = select((t0<=t) && (t<=t1), atan2(ell[6*1+4],ell[6*1+3]))$
  v=mod(v1-v2, 2*pi)$
  v =?(v<0):v+2*pi:v$
  v_npi = mod(v1-v2, 2*pi)$
  v_npi =?(v_npi>pi):v_npi-2*pi:v_npi$
  v_npi =?(v_npi<-pi):v_npi+2*pi:v_npi$
  j=4$
  vt[j][k] = t1$
  vmin[j][k]=min(v)$
  vmax[j][k]=max(v)$
  vmoy[j][k]=sum(v)/size(v)$
  j=5$
  vmin[j][k]=min(v_npi)$
  vmax[j][k]=max(v_npi)$
  vmoy[j][k]=sum(v_npi)/size(v_npi)$
 t0=t1$
};

 for j=4 to 5 {
    diffmin=max(abs(minmax_alp[3*(j-1)+1][q]-vmin[j]))$
    diffmax=max(abs(minmax_alp[3*(j-1)+3][q]-vmax[j]))$
    diffmoy=max(abs(minmax_alp[3*(j-1)+2][q]-vmoy[j]))$
    if ((diffmin>1E-12) || (diffmax>1E-12) || (diffmoy>1E-12)) then 
    {
        error("erreur min/max de pi_p(1)-pi_p(2)");
    };
 };


};
msg "min/max/moy a_p(1)-a_p(2), ... [-pi,pi] et [0,2*pi] .....................ok";
t=0,10$
write("extractok.dat", t);
}
catch {
msg "-----------------------------------------";
msg "-----------------------------------------";
msg "erreur dans la determination du min/max de a_p(1)-a_p(2) [-pi,pi] et [0,2*pi] ";
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


