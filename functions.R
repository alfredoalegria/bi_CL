####################################################
# Correlation functions employed in the manuscript
####################################################
corr <- function(choice,x,range){ 

    if(choice==1){corr = exp(-3*x/range);}           
    if(choice==2){corr = exp(-4.761905*x/range)*(1+4.761905*x/range);}
    if(choice==3){corr = 1/(1 + (4.358899*x/range)^2);}         
    return(corr);	
}
####################################################
# The 'partition' function returns two lists of indices, 'set1' and 'set2', 
# corresponding to labels 'a' and 'b' in the formulation of the method. The 
# inputs cx and cy corespond to the x and y coord. of the original grid, resp.
####################################################
partition <- function(m,cx,cy){

    set1 = c(); 
    set2 = c();    
    c = cbind(runif(m),runif(m));  # Generate nodes (blocks are formed around them)    
    ind = c();   # List containing the indices of the original grid already asigned to a block (initially empty)
    
    for(j in 1:m){
        cx[ind] = 3; # Once a point is asigned to a block, it cannot be selected again; to enforce this  
        cy[ind] = 3; # rule, we position the chosen point far from the domain, e.g., at coord. (3,3)
        
        aux = sqrt((cx-c[j,1])^2 + (cy-c[j,2])^2);
        ss = sort(aux);
  
        ind1 = which(aux==ss[1]); # The indices of the two closest points to the j-th node are identified
        ind2 = which(aux==ss[2]); 
  
        set1[j] = ind1;
        set2[j] = ind2;
  
        ind = c(c(ind),c(ind1,ind2)); # Systematically save the indices that have already been identified
    }
    return(cbind(set1,set2));  
}
####################################################
