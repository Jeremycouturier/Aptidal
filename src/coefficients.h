#ifndef _COEFFICIENTS_H_
#define _COEFFICIENTS_H_


/******** Coefficients of the Hamiltonian ********/
extern typ C_ppp1_1, C_ppp1_2, C_ppp1_3, C_ppp1_4, C_ppp1_5, C_ppp1_6, C_ppp1_7, C_ppp1_8, C_ppp1_9, C_ppp1_10, C_ppp1_11, C_ppp1_12, C_ppp1_13, C_ppp1_14, C_ppp1_15; //Resonance p:p+1
extern typ C_ppp2_1, C_ppp2_2, C_ppp2_3; //Resonance p:p+2
extern typ C_ppp3_1, C_ppp3_2, C_ppp3_3, C_ppp3_4; //Resonance p:p+3
extern typ C_00_1, C_00_2, C_00_3, C_00_4, C_00_5, C_00_6, C_00_7, C_00_8, C_00_9, C_00_10; //Resonance 0:0

/******** Array of pointers towards the functions resonance_pq ********/
extern void (*resonances[10][10])(typ alp, typ mi);



void resonance_init();


void resonance_00(typ alp, typ mi);


void resonance_12(typ alp, typ mi);


void resonance_23(typ alp, typ mi);


void resonance_34(typ alp, typ mi);


void resonance_45(typ alp, typ mi);


void resonance_56(typ alp, typ mi);


void resonance_67(typ alp, typ mi);


void resonance_78(typ alp, typ mi);


void resonance_89(typ alp, typ mi);


void resonance_13(typ alp, typ mi);


void resonance_35(typ alp, typ mi);


void resonance_57(typ alp, typ mi);


void resonance_79(typ alp, typ mi);


void resonance_14(typ alp, typ mi);


void resonance_25(typ alp, typ mi);


void resonance_47(typ alp, typ mi);


void resonance_58(typ alp, typ mi);


#endif

