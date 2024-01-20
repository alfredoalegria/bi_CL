#########################################
# Use the terminal to reach your directory. To compile the C code,
# enter 'R CMD SHLIB bi_cl.c'. This process will produce a shared 
# object file that can be dynamically linked to R
#########################################
source("functions.R");
source("bi_cl.R");
dyn.load("bi_cl.so");
library("geoR");
#########################################
# Generation of spatial sites
#########################################
set.seed(13);
sq = seq(0,1,by=0.03);
grid_sq = expand.grid(sq,sq); # Regular grid
num = nrow(grid_sq);
grid_sq = grid_sq+cbind(runif(num,-0.01,0.01), # Uniform perturbation
               runif(num,-0.01,0.01));
N = 500;
samp = sample(1:num,N);   # Sampling N locations
coordx = grid_sq[samp,1];
coordy = grid_sq[samp,2];
grid = cbind(coordx,coordy);
#########################################
# Simulation of a Gaussian process through Cholesky decomp.
#########################################
choice = 1;   # Model: Exponential=1; Matérn(1.5)=2; Cauchy=3
s2 = 1;      # Set variance 
range = 0.1; # Set range 
tau2 = 0.1;  # Set nugget 

dista = as.matrix(dist(grid,diag=TRUE,upper=TRUE));
covmat = s2*corr(choice,dista,range)+tau2*diag(N);
chol_cov = chol(covmat);
z = t(chol_cov)%*%rnorm(N);
#########################################
# Labeling half of the sites with superscript "a" and the rest
# with superscript "b" through the auxiliary function 'partition'.
# Build the bi-CL obj. funct. adding a certain number of config.
#########################################
m = N/2;
num_config = 5;  # Number of config.
set_a = array(dim=c(m,num_config)); 
set_b = array(dim=c(m,num_config)); 

for(r in 1:num_config){
   part = partition(m,coordx,coordy);
   set_a[,r] = part[,1]; 
   set_b[,r] = part[,2];
}

ds = 0.1;    # Fix threshold 'ds' for 0/1 weights
BICL <- function(par){
  aux = 0;
  for(r in 1:num_config){
    aux = aux + bi_cl(par,coordx[set_a[,r]],coordy[set_a[,r]],
                   coordx[set_b[,r]],coordy[set_b[,r]],
                   z[set_a[,r]],z[set_b[,r]],m,choice,ds);
  }	
  return(aux);
}
#########################################
# Estimating the parameters: define a starting point 'ini' 
# for the optimization algorithm; bi_cl uses an exponential  
# parameterization to guarantee positivity of the estimates
#########################################
ini = log(c(2,0.1,0.01)); 
EST <- optim(ini,BICL,control=list(reltol=1e-14,maxit=1e4,trace=TRUE));
print(exp(EST$par));

