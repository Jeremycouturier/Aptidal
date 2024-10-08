#programme qui merge les fichiers mpi
# necessite un dossier et le radical issu de intplastat??.par
# e.g.
# mergempi.sh DATA_saban  4planJN_tauscan
ls $1/$2.par || ( echo fichier inexistant : merge impossible && exit 1)
RAD=`cd $1 && ls $2_proc000.control | sed -e 's/_proc000//' | sed -e 's/.control$//' ` 
(cat $1/$2_proc*.control | sort -s -k 1,1 > $1/$RAD.control) || exit 1
(cat $1/$2_proc*.control_part | sort -s -k 1,1 > $1/$RAD.control_part) || exit 1
(cat $1/$2_proc*.ci_pla | sort -s -k 1,1 > $1/$RAD.ci_pla) || exit 1
(cat $1/$2_proc*.ci_part | sort -s -k 1,1 > $1/$RAD.ci_part) || exit 1
cat $1/$2_proc*.car_part | sort -s -k 1,1 > $1/$RAD.car_part
cat $1/$2_proc*.ell_part | sort -s -k 1,1 > $1/$RAD.ell_part
cat $1/$2_proc*.minmax_aei_part | sort -s -k 1,1 > $1/$RAD.minmax_aei_part
cat $1/$2_proc*.minmax_diffalp_part | sort -s -k 1,1 > $1/$RAD.minmax_diffalp_part
cat $1/$2_proc*.naf_alkhqp_part | sort -s -k 1,1 > $1/$RAD.naf_alkhqp_part
cat $1/$2_proc*.naf_alkh_part | sort -s -k 1,1 > $1/$RAD.naf_alkh_part
cat $1/$2_proc*.naf_diffalp_part | sort -s -k 1,1 > $1/$RAD.naf_diffalp_part
cat $1/$2_proc*.naf_alkhqp | sort -s -k 1,1 > $1/$RAD.naf_alkhqp
cat $1/$2_proc*.naf_alkh | sort -s -k 1,1 > $1/$RAD.naf_alkh
