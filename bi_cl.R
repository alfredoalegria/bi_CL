bi_cl <- function(par,cx_a,cy_a,cx_b,cy_b,z_a,z_b,m,choice,ds){
   
  storage.mode(cx_a) <- "double" 
  storage.mode(cy_a) <- "double" 
  storage.mode(cx_b) <- "double" 
  storage.mode(cy_b) <- "double" 
  storage.mode(z_a) <- "double"  
  storage.mode(z_b) <- "double"  
  storage.mode(par) <- "double" 
  
  o <- .C("bi_cl", choice=as.integer(choice),z_a=z_a,z_b=z_b,cx_a=cx_a,cy_a=cy_a,cx_b=cx_b,cy_b=cy_b,
             m=as.integer(m),par=par,ds=as.double(ds),sum=0)
  return(o$sum) 

}
