/* Evaluation of bi-CL */
#include <stdio.h>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>
#include "corr_models.h"

void bi_cl(int *choice,double *z_a,double *z_b,double *cx_a,double *cy_a,double *cx_b,double *cy_b,
       int *m,double *par,double *ds,double *sum)
{  
  double det1,det2,zi_a,zj_a,zi_b,zj_b,aux,s2,range,tau2,P,Q,R,S,m1,m2,W11,W12,W22,
       rij_aa,rij_bb,rij_ab,rij_ba,rii_ab,rjj_ab,
       d_aa_ij,d_bb_ij,d_ab_ij,d_ba_ij,d_ab_ii,d_ab_jj;
  
  s2 = exp(par[0]);   
  range = exp(par[1]);  
  tau2 = exp(par[2]); 
  
  for(int i=0;i<*m;i++){  
     for(int j=0;j<*m;j++){
		   
          d_aa_ij = hypot(cx_a[i]-cx_a[j],cy_a[i]-cy_a[j]);                

          if(i!=j & d_aa_ij<*ds){
                        
              zi_a = z_a[i]; zi_b = z_b[i];
              zj_a = z_a[j]; zj_b = z_b[j];
                          
              d_bb_ij = hypot(cx_b[i]-cx_b[j],cy_b[i]-cy_b[j]);  
              d_ab_ij = hypot(cx_a[i]-cx_b[j],cy_a[i]-cy_b[j]);  
              d_ba_ij = hypot(cx_b[i]-cx_a[j],cy_b[i]-cy_a[j]); 
              d_ab_ii = hypot(cx_a[i]-cx_b[i],cy_a[i]-cy_b[i]);  
              d_ab_jj = hypot(cx_a[j]-cx_b[j],cy_a[j]-cy_b[j]);           

              rij_aa = corr_model(*choice,d_aa_ij,range);  
              rij_bb = corr_model(*choice,d_bb_ij,range);  
              rij_ab = corr_model(*choice,d_ab_ij,range);  
              rij_ba = corr_model(*choice,d_ba_ij,range);  
              rii_ab = corr_model(*choice,d_ab_ii,range);  
              rjj_ab = corr_model(*choice,d_ab_jj,range);  
                                     
              det1 = pow(s2+tau2,2.0)-pow(s2*rjj_ab,2.0);                           
              P = s2*((s2+tau2)*rij_aa-s2*rij_ab*rjj_ab)/det1;
              Q = s2*((s2+tau2)*rij_ab-s2*rij_aa*rjj_ab)/det1;
              R = s2*((s2+tau2)*rij_ba-s2*rij_bb*rjj_ab)/det1;
              S = s2*((s2+tau2)*rij_bb-s2*rij_ba*rjj_ab)/det1;
                        
              m1 = P*zj_a+Q*zj_b;
              m2 = R*zj_a+S*zj_b;                      
                                                                                      
              W11 = tau2+s2*(1-P*rij_aa-Q*rij_ab);
              W12 = s2*(rii_ab-P*rij_ba-Q*rij_bb);     
              W22 = tau2+s2*(1-R*rij_ba-S*rij_bb);
                        
              det2 = W11*W22-pow(W12,2.0);
              aux = W22*pow(zi_a-m1,2.0)+W11*pow(zi_b-m2,2.0)-2*W12*(zi_a-m1)*(zi_b-m2);
                        
              *sum = *sum + 0.5*log(det2) + 0.5*aux/det2;   

         }
     }
  }
  
 
}







