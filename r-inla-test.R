library(INLA)

N = 100 #500, 5000, 25000, 100000  
x = rnorm(N, mean=6,sd=2)
y = rnorm(N, mean=x,sd=1) 
data = list(x=x,y=y,N=N) 

test <- inla(y~x, family = c("gaussian"), data = data, control.predictor=list(link=1)) 



