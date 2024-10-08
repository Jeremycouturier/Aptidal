/* programme pour generer le fichier tests_intpartdiststarn.par */
mu=0.2958122082855911D-03;
ell=vnumR[30,0,0.6,0,0,0.2];
t=0,0.8,0.05;
for j=1 to 6 { ell[j]=ell[j]+t*0$ };
ell[3]=t;
ell[1][:]=30;
ell[5][:]=2E-2*cos(t);
ell[6][:]=2E-2*sin(t);
ell_to_xyz(ell,mu,x,xp);
w=1,size(t);
write(tests_filepartcidistn.dat,"PART%g 5 "+6*"%g "+"\n", w, x, xp);
!"cat tests_filepartcidistn.dat";



/* programme pour generer le fichier tests_intpartdistplan.par */
mu=0.2958122082855911D-03;
ell=vnumR[30,0,0.6,0,0,0.2];
t=0.05,0.1,0.003;
for j=1 to 6 { ell[j]=ell[j]+t*0$ };
ell[3]=t;
ell[1][:]=30;
ell[5][:]=2E-2*cos(t);
ell[6][:]=2E-2*sin(t);
ell_to_xyz(ell,mu,x,xp);
w=1,size(t);
write(tests_filepartcidistplan.dat,"PART%g 5 "+6*"%g "+"\n", w, x, xp);
!"cat tests_filepartcidistplan.dat";
