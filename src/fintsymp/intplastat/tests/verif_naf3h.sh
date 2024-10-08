#!/bin/bash
rm -r -f DATA_test2 DATA_test
cp tests_fileci.dat tests_fileci1.dat
mkdir DATA_test
../intplastat.x  tests_intplastatnaf3h.par >/dev/null || exit 1

#verification des analyses en frequence
rm -f extractok.dat
(trip<<__EOF
try {
npla=4$
nfreq=5$
fnaf_alkhqp="./DATA_test/test_proc000.naf_diffalp"$
fmt="%s %g "+6*nfreq*" %g"$
vnumR t_naf_alkhqp,naf_alkhqp[1:6*nfreq]$
read(fnaf_alkhqp, fmt, rad_naf_alkhqp, t_naf_alkhqp, naf_alkhqp);
_info off$
for ci=1 to 3 {
ci_cur=ci$
rad=str("N000%d", ci)$
fell="./DATA_test/test_"+rad+".ell"$
vnumR t, ell[1:npla*6];
read(fell, t, ell);

id=1,size(t_naf_alkhqp)$
q=dimtovnumR(select(rad_naf_alkhqp==rad,vnumtodim(id)))$

tsel=vnumR[0:18E1:36E1]$


vnumR vfsrc([1:5])$
vnumC vzsrc([1:5])$

kpla=0$
t0=tsel[1]$
kpla_cur = kpla$
for k=1 to size(tsel)-1
{
 //msg("kpla=%d k=%d\n", kpla, k);
 t1 = tsel[k+1]$
 vt=select((t0<=t) && (t<=t1), t)$
 vr = select((t0<=t) && (t<=t1), cos(ell[6*(kpla)+2]-ell[6*(kpla+1)+2]))$
 vi = select((t0<=t) && (t<=t1), sin(ell[6*(kpla)+2]-ell[6*(kpla+1)+2]))$
  naftab2(vr,vi,A,F,size(vr),t[2]-t[1], t0, nfreq);
  for j=1 to nfreq { vfsrc[j] = naf_alkhqp[kpla*6*nfreq+(1-1)*3*nfreq+3*(j-1)+1][q[k]]$ };
  for j=1 to nfreq { vzsrc[j] = naf_alkhqp[kpla*6*nfreq+(1-1)*3*nfreq+3*(j-1)+2][q[k]]+I*naf_alkhqp[kpla*6*nfreq+(1-1)*3*nfreq+3*(j-1)+3][q[k]]$ };
  /*writes(F-vfsrc,F,vfsrc);
  writes(vzsrc-A, A, vzsrc);
  S=sertrig(A,F,x);
  plot(vt, vr-real(evalnum(S,COMPLEX,(x,vt))));
  S=sertrig(vzsrc,vfsrc,x);
  replot(vt, vr-real(evalnum(S,COMPLEX,(x,vt))));*/
  difffreq=max(abs(F-vfsrc))$
  diffzamp=max(abs(vzsrc-A))$
  if (difffreq>5E-9) then {  error("erreur freq de a*expi(l)"); };
  if (diffzamp>5E-9) then {  error("erreur zamp de a*expi(l)"); };

  t1 = tsel[k+1]$
  vt=select((t0<=t) && (t<=t1), t)$
  pi1=atan2(ell[6*kpla+4],ell[6*kpla+3])$
  pi2=atan2(ell[6*(kpla+1)+4],ell[6*(kpla+1)+3])$
  vr = select((t0<=t) && (t<=t1), cos(pi1-pi2))$
  vi = select((t0<=t) && (t<=t1), sin(pi1-pi2))$
  naftab2(vr,vi,A,F,size(vr),t[2]-t[1], t0, nfreq);
  for j=1 to nfreq { vfsrc[j] = naf_alkhqp[kpla*6*nfreq+(2-1)*3*nfreq+3*(j-1)+1][q[k]]$ };
  for j=1 to nfreq { vzsrc[j] = naf_alkhqp[kpla*6*nfreq+(2-1)*3*nfreq+3*(j-1)+2][q[k]]+I*naf_alkhqp[kpla*6*nfreq+(2-1)*3*nfreq+3*(j-1)+3][q[k]]$ };
  
  /*writes(F-vfsrc,F,vfsrc);
  writes(vzsrc-A, A, vzsrc);
  S=sertrig(A,F,x);
  plot(vt, vr-real(evalnum(S,COMPLEX,(x,vt))));
  S=sertrig(vzsrc,vfsrc,x);
  replot(vt, vr-real(evalnum(S,COMPLEX,(x,vt))));*/
  difffreq=max(abs(F-vfsrc))$
  diffzamp=max(abs(vzsrc-A))$
  if (difffreq>5E-9) then {  error("erreur freq de k+i*h"); };
  if (diffzamp>5E-9) then {  error("erreur zamp de k+i*h"); };

 t0=t1$
 };
};
msg "naf exp(i*(l1-l2)),.....................ok";
t=0,10$
write("extractok.dat", t);
}
catch {
msg "-----------------------------------------";
msg "-----------------------------------------";
msg "erreur dans la determination du naf exp(i*(l1-l2))";
ci_cur;
difffreq;
diffzamp;
writes(F-vfsrc,F,vfsrc);
writes(vzsrc-A, A, vzsrc);
S=sertrig(A,F,x);
/*  plot(vt, vr);
  replot(vt, vi);
  pause(10);*/
  
  S=sertrig(vzsrc,vfsrc,x);
  plot(vt, vr-real(evalnum(S,COMPLEX,(x,vt))));
  replot(vt, vi-imag(evalnum(S,COMPLEX,(x,vt))));
  S=sertrig(vzsrc,vfsrc,x);
  replot(vt, vr-real(evalnum(S,COMPLEX,(x,vt))));
  replot(vt, vi-imag(evalnum(S,COMPLEX,(x,vt))));
  pause(20);
msg "-----------------------------------------";
error("-----------------------------------------");
};
__EOF
) | grep -v BIENVENUE | grep -v Taper | grep -v apprecie | sed '/^$/d'
    ls extractok.dat >/dev/null ||  exit 1 
    rm -f extractok.dat
rm -r -f DATA_test2 DATA_test
exit 0


