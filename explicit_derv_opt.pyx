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
rcm1 = N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31));
rcm2 = N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32));
rcm3 = N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33));
rcm1_sq = rcm1 * rcm1;
rcm2_sq = rcm2 * rcm2;
rcm3_sq = rcm3 * rcm3;
rcm_sq = rcm1_sq + rcm2_sq + rcm3_sq;
rcm = sqrt(rcm_sq);
t34 = -rcm3 / rcm / rcm_sq * rcm1 * resid_N;

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

cosb_N2 = -rcm3 / rcm / rcm_sq * rcm2 * resid_N;

## Derive cosb wrt N3
t3 = 1 / (m_N + 3 * m_H);
t5 = -t3 * m_N + 1;
t12 = pow((N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31))), 2.0);
t19 = pow((N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32))), 2.0);
t26 = pow((N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33))), 2.0);
t27 = t12 + t19 + t26;
t28 = sqrt(t27);
t35 = t5 / t28 - t26 / t28 / t27 * t5;

cosb_N3 = resid_N / rcm - rcm3_sq / rcm / rcm_sq * resid_N;

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

cosb_H11 = rcm3 / rcm / rcm_sq * rcm1 * resid_H;

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

cosb_H12 = rcm3 / rcm / rcm_sq * rcm2 * resid_H;

## Derive cosb wrt H13
t3 = 1 / (m_N + 3 * m_H);
t4 = t3 * m_H;
t11 = pow((N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31))), 2);
t18 = pow((N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32))), 2);
t25 = pow((N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33))), 2);
t26 = t11 + t18 + t25;
t27 = sqrt(t26);
t34 = -t4 / t27 + t25 / t27 / t26 * t4;

cosb_H13 = -resid_H / rcm + rcm3_sq / rcm / rcm_sq * resid_H;

## Derive cosb wrt H21
t3 = 1 / (m_N + 3 * m_H);
t9 = N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33));
t15 = N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31));
t16 = t15 * t15;
t23 = (int) pow((double) (N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32))), (double) 2);
t24 = t9 * t9;
t25 = t16 + t23 + t24;
t26 = sqrt((double) t25);
t32 = (double) t9 / t26 / (double) t25 * (double) t15 * (double) t3 * (double) m_H;

cosb_H21 =  rcm3 / rcm / rcm_sq * rcm1 * resid_H;

## Derive cosb wrt H22
t3 = 1 / (m_N + 3 * m_H);
t9 = N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33));
t16 = (int) pow((double) (N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31))), (double) 2);
t22 = N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32));
t23 = t22 * t22;
t24 = t9 * t9;
t25 = t16 + t23 + t24;
t26 = sqrt((double) t25);
t32 = (double) t9 / t26 / (double) t25 * (double) t3 * (double) t22 * (double) m_H;

cosb_H22 = rcm3 / rcm / rcm_sq * rcm2 * resid_H;

## Derive cosb wrt H23
t3 = 1 / (m_N + 3 * m_H);
t4 = t3 * m_H;
t11 = (int) pow((double) (N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31))), (double) 2);
t18 = (int) pow((double) (N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32))), (double) 2);
t25 = (int) pow((double) (N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33))), (double) 2);
t26 = t11 + t18 + t25;
t27 = sqrt((double) t26);
t34 = -(double) t4 / t27 + (double) t25 / t27 / (double) t26 * (double) t4;

cosb_H23 = -resid_H / rcm  + rcm3_sq / rcm / rcm_sq * resid_H;

## Derive cosb wrt H31
t3 = 1 / (m_N + 3 * m_H);
t9 = N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33));
t15 = N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31));
t16 = t15 * t15;
t23 = (int) pow((double) (N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32))), (double) 2);
t24 = t9 * t9;
t25 = t16 + t23 + t24;
t26 = sqrt((double) t25);
t32 = (double) t9 / t26 / (double) t25 * (double) t3 * (double) t15 * (double) m_H;

cosb_H31 = rcm3 / rcm / rcm_sq * rcm1 * resid_H;

## Derive cosb wrt H32
t3 = 1 / (m_N + 3 * m_H);
t9 = N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33));
t16 = (int) pow((double) (N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31))), (double) 2);
t22 = N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32));
t23 = t22 * t22;
t24 = t9 * t9;
t25 = t16 + t23 + t24;
t26 = sqrt((double) t25);
t32 = (double) t9 / t26 / (double) t25 * (double) t3 * (double) t22 * (double) m_H;

cosb_H32 = rcm3 / rcm / rcm_sq * rcm2 * resid_H;

## Derive cosb wrt H33
t3 = 1 / (m_N + 3 * m_H);
t4 = t3 * m_H;
t11 = (int) pow((double) (N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31))), (double) 2);
t18 = (int) pow((double) (N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32))), (double) 2);
t25 = (int) pow((double) (N3 - t3 * (m_N * N3 + m_H * (H13 + H23 + H33))), (double) 2);
t26 = t11 + t18 + t25;
t27 = sqrt((double) t26);
t34 = -(double) t4 / t27 + (double) t25 / t27 / (double) t26 * (double) t4;

cosb_H33 = -resid_H / rcm + rcm3_sq / rcm / rcm_sq * resid_H;

## Derive a wrt N1
t3 = 1 / (m_N + 3 * m_H);
t5 = -t3 * m_N + 1;
t12 = (int) pow((double) (N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31))), (double) 2);
t19 = (int) pow((double) (N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32))), (double) 2);
t20 = t12 + t19;
t21 = sqrt((double) t20);
t32 = sqrt((double) (1 - t12 / t20));
t35 = -((double) t5 / t21 - (double) t12 / t21 / (double) t20 * (double) t5) / t32;

a_N1 =  -(resid_N / rcm12 - rcm1_sq / rcm12 / rcm12_sq * resid_N) / (sqrt(1 - rcm1_sq / rcm12_sq));

## Derive a wrt N2
t3 = 1 / (m_N + 3 * m_H);
t9 = N1 - t3 * (m_N * N1 + m_H * (H11 + H21 + H31));
t10 = t9 * t9;
t16 = N2 - t3 * (m_N * N2 + m_H * (H12 + H22 + H32));
t17 = t16 * t16;
t18 = t10 + t17;
t19 = sqrt((double) t18);
t29 = sqrt((double) (1 - t10 / t18));
t32 = (double) t9 / t19 / (double) t18 * (double) t16 * (double) (-t3 * m_N + 1) / t29;

a_N2 = rcm1 / rcm12 / rcm12_sq * rcm2 * resid_N / a_arccos_deriv;