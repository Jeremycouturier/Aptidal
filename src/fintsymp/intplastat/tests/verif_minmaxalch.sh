#!/bin/bash
rm -r -f DATA_test2 DATA_test
cp tests_fileci.dat tests_fileci1.dat
mkdir DATA_test
../intplastat.x  tests_intplastatminmaxalch.par >/dev/null || exit 1

#verification des min/max a_p(1)-a_p(2), ... avec angle continu
rm -f extractok.dat
(trip<<__EOF
try {
npla=2$
fminmax_alc="./DATA_test/test_proc000.minmax_diffalc"$
fmt="%s %g "+9*" %g"$
vnumR t_minmax_alc,minmax_alc[1:9]$
read(fminmax_alc, fmt, rad_minmax_alc, t_minmax_alc, minmax_alc);
for ci=1 to 3 {
rad=str("N000%d", ci)$
fell="./DATA_test/test_"+rad+".ell"$
vnumR t, ell[1:npla*6];
read(fell, t, ell);

id=1,size(t_minmax_alc)$
q=dimtovnumR(select(rad_minmax_alc==rad,vnumtodim(id)))$

macro angle[_TANG]{
/*********************************************************/
/* macro qui retrouve la continuite d'un  tableau d'angles */
/*********************************************************/
  private _n, _Dang, _saut, _csaut;
  _n=size(_TANG)$
  _Dang = _TANG[2:_n]-_TANG[1:_n-1]$
  _saut = ?(_Dang>2):-1:(?(_Dang<-2):1:0)$
  _csaut = accum(_saut)$
   _TANG[2:] = _TANG[2:] + _csaut*2*pi$
};

vnumR vt[1:5]([1:size(q)])$
vnumR vmin[1:5]([1:size(q)])$
vnumR vmax[1:5]([1:size(q)])$
vnumR vmoy[1:5]([1:size(q)])$
t0=0$
for k=1 to size(q)
{
 t1 = t_minmax_alc[q[k]]$
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
    diffmin=max(abs(minmax_alc[3*(j-1)+1][q]-vmin[j]))$
    diffmax=max(abs(minmax_alc[3*(j-1)+3][q]-vmax[j]))$
    diffmoy=max(abs(minmax_alc[3*(j-1)+2][q]-vmoy[j]))$
    if ((diffmin>1E-12) || (diffmax>1E-12) || (diffmoy>1E-12)) then 
    {
        error("erreur min/max de a_p(1)-a_p(2) avec angle continu ");
    };
 };

t0=0$
for k=1 to size(q)
{
  t1 = t_minmax_alc[q[k]]$
  v1 = select((t0<=t) && (t<=t1), ell[6*0+2])$
  v2 = select((t0<=t) && (t<=t1), ell[6*1+2])$
  v=v1-v2$
  %angle[[v]];
  j=2$
  vt[j][k] = t1$
  vmin[j][k]=min(v)$
  vmax[j][k]=max(v)$
  vmoy[j][k]=sum(v)/size(v)$
  t0=t1$
};

 for j=2 to 2 {
    diffmin=max(abs(minmax_alc[3*(j-1)+1][q]-vmin[j]))$
    diffmax=max(abs(minmax_alc[3*(j-1)+3][q]-vmax[j]))$
    diffmoy=max(abs(minmax_alc[3*(j-1)+2][q]-vmoy[j]))$
    if ((diffmin>1E-12) || (diffmax>1E-12) || (diffmoy>1E-12)) then 
    {
        error("erreur min/max de l_p(1)-l_p(2) avec angle continu ");
    };
 };


t0=0$
for k=1 to size(q)
{
  t1 = t_minmax_alc[q[k]]$
  v1 = select((t0<=t) && (t<=t1), atan2(ell[6*0+4],ell[6*0+3]))$
  v2 = select((t0<=t) && (t<=t1), atan2(ell[6*1+4],ell[6*1+3]))$
  v=v1-v2$
  %angle[[v]];
  j=3$
  vt[j][k] = t1$
  vmin[j][k]=min(v)$
  vmax[j][k]=max(v)$
  vmoy[j][k]=sum(v)/size(v)$
 t0=t1$
};

 for j=3 to 3 {
    diffmin=max(abs(minmax_alc[3*(j-1)+1][q]-vmin[j]))$
    diffmax=max(abs(minmax_alc[3*(j-1)+3][q]-vmax[j]))$
    diffmoy=max(abs(minmax_alc[3*(j-1)+2][q]-vmoy[j]))$
    if ((diffmin>1E-12) || (diffmax>1E-12) || (diffmoy>1E-12)) then 
    {
        error("erreur min/max de pi_p(1)-pi_p(2) avec angle continu ");
    };
 };


};
msg "min/max/moy a_p(1)-a_p(2),...  avec angle continu.....................ok";
t=0,10$
write("extractok.dat", t);
}
catch {
msg "-----------------------------------------";
msg "-----------------------------------------";
msg "erreur dans la determination du min/max de a_p(1)-a_p(2) avec angle continu ";
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


