#cdef double cosb_n1(double )
## REFERENCE Derive cosb wrt N1
t3 = 1 / (m_N + 3 * m_H);
t9 = N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33));
t15 = N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31));
t16 = t15 * t15;
t23 = pow (N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32))), 2.0);
t24 = t9 * t9;
t25 = t16 + t23 + t24;
t26 = sqrt t25);
t34 =  t9 / t26 / t25 * t15 * resid_N;

## Derive cosb wrt N1
t3 = 1 / (m_N + 3 * m_H);
resid_N = -t3 * m_N + 1;
Ncm1 = N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31));
Ncm2 = N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32));
Ncm3 = N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33));
Ncm1_sq = Ncm1 * Ncm1;
Ncm2_sq = Ncm2 * Ncm2;
Ncm3_sq = Ncm3 * Ncm3;
Ncm_sq = Ncm1_sq + Ncm2_sq + Ncm3_sq;
Ncm = sqrt(Ncm_sq);
t34 = -Ncm3 / Ncm / Ncm_sq * Ncm1 * resid_N;

## Derive cosb wrt N2
t3 = 1 / (m_N + 3 * m_H);
t9 = N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33));
t16 = pow((N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31))), 2.0);
t22 = N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32));
t23 = t22 * t22;
t24 = t9 * t9;
t25 = t16 + t23 + t24;
t26 = sqrt(t25);
t34 = -t9 / t26 / t25 * t22 * resid_N;

cosb_N2 = -Ncm3 / Ncm / Ncm_sq * Ncm2 * resid_N;

## Derive cosb wrt N3
t3 = 1 / (m_N + 3 * m_H);
t5 = -t3 * m_N + 1;
t12 = pow((N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31))), 2.0);
t19 = pow((N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32))), 2.0);
t26 = pow((N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33))), 2.0);
t27 = t12 + t19 + t26;
t28 = sqrt(t27);
t35 = t5 / t28 - t26 / t28 / t27 * t5;

cosb_N3 = resid_N / Ncm - Ncm3_sq / Ncm / Ncm_sq * resid_N;

## Derive cosb wrt H11
t3 = 1 / (m_N + 3 * m_H);
resid_H = t3 * m_H	
t9 = N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33));
t15 = N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31));
t16 = t15 * t15;
t23 = pow (N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32))), 2);
t24 = t9 * t9;
t25 = t16 + t23 + t24;
t26 = sqrt(t25);
t32 = t9 / t26 / t25 * t15 * t3 * m_H;

cosb_H11 = Ncm3 / Ncm / Ncm_sq * Ncm1 * resid_H;

## Derive cosb wrt H12
t3 = 1 / (m_N + 3 * m_H);
t9 = N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33));
t16 = pow((N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31))), 2);
t22 = N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32));
t23 = t22 * t22;
t24 = t9 * t9;
t25 = t16 + t23 + t24;
t26 = sqrt(t25);
t32 = t9 / t26 / t25 * t22 * m_H * t3;

cosb_H12 = Ncm3 / Ncm / Ncm_sq * Ncm2 * resid_H;

## Derive cosb wrt H13
t3 = 1 / (m_N + 3 * m_H);
t4 = t3 * m_H;
t11 = pow((N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31))), 2);
t18 = pow((N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32))), 2);
t25 = pow((N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33))), 2);
t26 = t11 + t18 + t25;
t27 = sqrt(t26);
t34 = -t4 / t27 + t25 / t27 / t26 * t4;

cosb_H13 = -resid_H / Ncm + Ncm3_sq / Ncm / Ncm_sq * resid_H;

## Derive cosb wrt H21
t3 = 1 / (m_N + 3 * m_H);
t9 = N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33));
t15 = N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31));
t16 = t15 * t15;
t23 = (int) pow((N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32))), 2);
t24 = t9 * t9;
t25 = t16 + t23 + t24;
t26 = sqrt(t25);
t32 = t9 / t26 / t25 * t15 * t3 * m_H;

cosb_H21 =  Ncm3 / Ncm / Ncm_sq * Ncm1 * resid_H;

## Derive cosb wrt H22
t3 = 1 / (m_N + 3 * m_H);
t9 = N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33));
t16 = (int) pow((N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31))), 2);
t22 = N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32));
t23 = t22 * t22;
t24 = t9 * t9;
t25 = t16 + t23 + t24;
t26 = sqrt(t25);
t32 = t9 / t26 / t25 * t3 * t22 * m_H;

cosb_H22 = Ncm3 / Ncm / Ncm_sq * Ncm2 * resid_H;

## Derive cosb wrt H23
t3 = 1 / (m_N + 3 * m_H);
t4 = t3 * m_H;
t11 = (int) pow((N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31))), 2);
t18 = (int) pow((N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32))), 2);
t25 = (int) pow((N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33))), 2);
t26 = t11 + t18 + t25;
t27 = sqrt(t26);
t34 = -t4 / t27 + t25 / t27 / t26 * t4;

cosb_H23 = -resid_H / Ncm  + Ncm3_sq / Ncm / Ncm_sq * resid_H;

## Derive cosb wrt H31
t3 = 1 / (m_N + 3 * m_H);
t9 = N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33));
t15 = N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31));
t16 = t15 * t15;
t23 = (int) pow((N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32))), 2);
t24 = t9 * t9;
t25 = t16 + t23 + t24;
t26 = sqrt(t25);
t32 = t9 / t26 / t25 * t3 * t15 * m_H;

cosb_H31 = Ncm3 / Ncm / Ncm_sq * Ncm1 * resid_H;

## Derive cosb wrt H32
t3 = 1 / (m_N + 3 * m_H);
t9 = N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33));
t16 = (int) pow((N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31))), 2);
t22 = N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32));
t23 = t22 * t22;
t24 = t9 * t9;
t25 = t16 + t23 + t24;
t26 = sqrt(t25);
t32 = t9 / t26 / t25 * t3 * t22 * m_H;

cosb_H32 = Ncm3 / Ncm / Ncm_sq * Ncm2 * resid_H;

## Derive cosb wrt H33
t3 = 1 / (m_N + 3 * m_H);
t4 = t3 * m_H;
t11 = (int) pow((N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31))), 2);
t18 = (int) pow((N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32))), 2);
t25 = (int) pow((N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33))), 2);
t26 = t11 + t18 + t25;
t27 = sqrt(t26);
t34 = -t4 / t27 + t25 / t27 / t26 * t4;

cosb_H33 = -resid_H / Ncm + Ncm3_sq / Ncm / Ncm_sq * resid_H;

## Derive a wrt N1
t3 = 1 / (m_N + 3 * m_H);
t5 = -t3 * m_N + 1;
t12 = (int) pow((N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31))), 2);
t19 = (int) pow((N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32))), 2);
t20 = t12 + t19;
t21 = sqrt(t20);
t32 = sqrt((1 - t12 / t20));
t35 = -(t5 / t21 - t12 / t21 / t20 * t5) / t32;

a_N1 =  -(resid_N / Ncm12 - Ncm1_sq / Ncm12 / Ncm12_sq * resid_N) / (sqrt(1 - Ncm1_sq / Ncm12_sq));

## Derive a wrt N2
t3 = 1 / (m_N + 3 * m_H);
t9 = N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31));
t10 = t9 * t9;
t16 = N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32));
t17 = t16 * t16;
t18 = t10 + t17;
t19 = sqrt(t18);
t29 = sqrt((1 - t10 / t18));
t32 = t9 / t19 / t18 * t16 * (-t3 * m_N + 1) / t29;

a_N2 = Ncm1 / Ncm12 / Ncm12_sq * Ncm2 * resid_N / a_arccos_deriv;

## Derive a wrt N3
a_N3 = 0.0

## Derive a wrt H11
t3 = 1 / (m_N + 3 * m_H);
t4 = t3 * m_H;
t11 = (int) pow((N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31))), 2);
t18 = (int) pow((N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32))), 2);
t19 = t11 + t18;
t20 = sqrt(t19);
t31 = sqrt((1 - t11 / t19));
t34 = -(-t4 / t20 + t11 / t20 / t19 * t4) / t31;

a_H11 = (resid_H / Ncm12 - Ncm1_sq / Ncm12 / Ncm12_sq * resid_H) / a_arccos_deriv;

## Deriev g wrt N1
t3 = 1 / (m_N + 3 * m_H);
t8 = t3 * (m_N * N2 + m_H * (H12 + H22 + H32));
t9 = N2 - t8;
t14 = t3 * (m_N * N3 + m_H * (H13 + H23 + H33));
t15 = H13 - t14;
t17 = N3 - t14;
t18 = H12 - t8;
t20 = t9 * t15 - t17 * t18;
t25 = t3 * (m_N * N1 + m_H * (H11 + H21 + H31));
t26 = N1 - t25;
t27 = t26 * t26;
t28 = t9 * t9;
t29 = t27 + t28;
t30 = 1 / t29;
t33 = sqrt((-t27 * t30 + 1));
t38 = -t3 * m_N + 1;
t41 = t29 * t29;
t51 = -t17 * t3 * m_N - t38 * t15;
t53 = sqrt(t29);
t54 = 0.1e1 / t53;
t57 = H11 - t25;
t59 = -t26 * t15 + t17 * t57;
t68 = t20 * t20;
t69 = t59 * t59;
t72 = t26 * t18 - t9 * t57;
t73 = t72 * t72;
t74 = t68 + t69 + t73;
t75 = sqrt(t74);
t81 = t59 * t26 * t54 - t20 * t33;
t95 = t81 * t81;
t99 = sqrt(0.1e1 - t95 / t74);
t102 = -((-t20 / t33 * (-t26 * t30 * t38 + t27 * t26 / t41 * t38) + t51 * t26 * t54 + t59 * t38 * t54 - t59 * t27 / t53 / t29 * t38) / t75 - t81 / t75 / t74 * (t59 * t51 + t72 * (t9 * t3 * m_N + t38 * t18))) / t99;


## Derive g wrt H33
t3 = 1 / (m_N + 3 * m_H);
t8 = t3 * (m_N * N2 + m_H * (H12 + H22 + H32));
t9 = N2 - t8; ################################################### Ncm2
t12 = t3 * m_H; ################################################# resid_H
t13 = H12 - t8; ################################################# H1cm2
t15 = -t9 * t3 * m_H + t12 * t13; #################################### comp1
t20 = t3 * (m_N * N1 + m_H * (H11 + H21 + H31));
t21 = N1 - t20; ################################################# Ncm1
t22 = t21 * t21; ################################################ Ncm1_sq
t23 = t9 * t9; ################################################## Ncm2_sq
t24 = t22 + t23; ################################################ Ncm12_sq
t28 = sqrt((1 - t22 / t24)); ########################### a_arccos_deriv
t32 = H11 - t20; ################################################ H1cm1
t34 = t21 * t3 * m_H - t12 * t32; #################################### comp2
t36 = sqrt(t24); ####################################### Ncm12
t37 = 0.1e1 / t36;
t44 = t3 * (m_N * N3 + m_H * (H13 + H23 + H33));
t45 = H13 - t44;
t47 = N3 - t44;
t49 = -t47 * t13 + t9 * t45;
t50 = t49 * t49;
t53 = -t21 * t45 + t47 * t32;
t54 = t53 * t53;
t58 = (int) pow((t21 * t13 - t9 * t32), 2);
t59 = t50 + t54 + t58;
t60 = sqrt(t59);
t66 = t53 * t21 * t37 - t49 * t28;
t76 = t66 * t66;
t80 = sqrt(0.1e1 - t76 / t59);
t83 = -((t34 * t21 * t37 - t15 * t28) / t60 - t66 / t60 / t59 * (t49 * t15 + t53 * t34)) / t80;

