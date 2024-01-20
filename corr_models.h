#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<R.h>
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <Rinternals.h>

double corr_model(int choice,double dista,double range);
double corr_model(int choice,double dista,double range) 
{
    double corr;
      
    // Exponential  
    if(choice==1){ 
		corr = exp(-3*dista/range);
	}             
    // Matérn (1.5)  
    if(choice==2){ 
		corr = exp(-4.761905*dista/range)*(1+4.761905*dista/range);
	}  	
    // Cauchy
    if(choice==3){ 
		corr = 1/(1 + pow(4.358899*dista/range,2.0));
	}   
	
    return(corr);
}























