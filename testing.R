library(Rcpp)
sourceCpp("testing.cpp")

print("C++ Function:")
system.time(dBvZINB4_Expt(x,y,a0,a1,a2,b1,b2,p1,p2,p3,p4)) 

fib_r <- function(n){
  if(n==1||n==2) return(1)
  return(fib_r(n-1)+fib_r(n-2))
}

print("R Function:")
system.time(fib_r(30))

