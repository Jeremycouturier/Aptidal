#programme qui merge les fichiers mpi
# necessite un dossier et le radical issu de intplastat??.par
# e.g.
# mergempi.sh DATA_saban  4planJN_tauscan
ls $1/$2.par || ( echo fichier iniexistant : merge impossible && exit 1)
RAD=`cd $1 && ls $2_proc000.control | sed -e 's/_proc000//' | sed -e 's/.control$//' ` 
(cat $1/$2_proc*.control | sort -s -k 1,1 > $1/$RAD.control) || exit 1
(cat $1/$2_proc*.ci | sort -s -k 1,1 > $1/$RAD.ci) || exit 1
cat $1/$2_proc*.minmax_aei | sort -s -k 1,1 > $1/$RAD.minmax_aei
cat $1/$2_proc*.minmax_diffalp | sort -s -k 1,1 > $1/$RAD.minmax_diffalp
cat $1/$2_proc*.naf_alkhqp | sort -s -k 1,1 > $1/$RAD.naf_alkhqp
cat $1/$2_proc*.naf_alkh | sort -s -k 1,1 > $1/$RAD.naf_alkh
cat $1/$2_proc*.naf_diffalp | sort -s -k 1,1 > $1/$RAD.naf_diffalp
