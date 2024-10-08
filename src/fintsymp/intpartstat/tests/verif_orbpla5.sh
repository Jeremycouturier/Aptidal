#!/bin/bash
echo "verification de if_orb_pla==1 .... en cours"

rm -r -f DATA_testnref DATA_test2
mkdir DATA_testnref
../intpartstat.x tests_intpartstatorbpla0.par >/dev/null || exit 1
mkdir DATA_test2
../intpartstat.x tests_intpartstatorbpla1.par >/dev/null || exit 1

rm -f extractok.dat
(trip <<__EOF
try { 
//particule
filename=["test_proc000.car_part": "test_proc000.ell_part"]$
ncomp=[6:6]$
for k=1 to size(filename)
{
    n=ncomp[k]$
    vnumR ts2all, t2all, xyz2all[1:n];
    read("DATA_test2/"+filename[k],ts2all, t2all,xyz2all);

    vnumR tsrall, trall, xyzrall[1:n];
    read("DATA_testnref/"+filename[k], tsrall, trall, xyzrall);

    for vipart=0 to 10 {
 
        t2=select(ts2all==vipart, t2all)$
        vnumR xyz2[1:n];
        for j=1 to n { xyz2[j]=select(ts2all==vipart, xyz2all[j])$ };
        vnumR xyzr[1:n];
        tr=select(tsrall==vipart, trall)$ tr=tr[::10]$
        for j=1 to n { xyzr[j]=select(tsrall==vipart, xyzrall[j])$ xyzr[j]=xyzr[j][::10]$ };

        for j=1 to n {
            if (max(abs(xyz2[j]-xyzr[j]))>=1E-13) then {
            j;
            max(abs(xyz2[j]-xyzr[j]));
            msg "-----------------------------------------";
            msg "-----------------------------------------";
            msg "erreur dans les coordonnees "+filename[k]+"  de la particule "+str(vipart);
            msg "-----------------------------------------";
            error("-----------------------------------------");
            };
        };
    };
};

//planete
filename=["test_proc000.ell" ]$
ncomp=[18]$
for k=1 to size(filename)
{
    n=ncomp[k]$
    vnumR t2all, xyz2all[1:n];
    read("DATA_test2/"+filename[k], t2all,xyz2all);

    vnumR trall, xyzrall[1:n];
    read("DATA_testnref/"+filename[k], trall, xyzrall);

    trall=trall[::10]$
    for j=1 to n { xyzrall[j]=xyzrall[j][::10]$ };

    for j=1 to n {
       if (max(abs(xyz2all[j]-xyzrall[j]))>=1E-11) then {
       j;
       max(abs(xyz2all[j]-xyzrall[j]));
       msg "-----------------------------------------";
       msg "-----------------------------------------";
       msg "erreur dans les coordonnees des planetes "+filename[k];
       msg "-----------------------------------------";
       error("-----------------------------------------");
       };
    };
};

t=0,10$
write("extract.dat", t);
} catch { msg "une erreur a eu lieu!!!";};
__EOF
) | grep -v BIENVENUE | grep -v Taper | grep -v apprecie | sed '/^$/d'
ls extract.dat >/dev/null || { echo "verification de if_orb_pla==1 .................... FAILED!" && exit 1 ; }
rm -f extract.dat

rm -r -f DATA_testnref DATA_test2
echo "verification de if_orb_pla==1 .................... ok"
exit 0


