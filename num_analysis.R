# Numerical Analysis 
#===================
library(numDeriv)
library(alabama)
library(readxl)
library(deSolve)
library(linprog)
library(ggplot2)
library(quadprog)
library(minpack.lm)
library(BB)
library(Rsolnp)
library(dplyr)
library(lpSolve)
library(gaussquad)
library(cubature)
library(ucminf)
library(tidyverse)
library(pracma)
options(digits=5)

#Working directory
script_name <- 'num_analysis.R'
ruta <- gsub(rstudioapi::getActiveDocumentContext()$path,pattern = script_name,replacement = '')
setwd(ruta)

# 1.- Numerical Differentiation
#------------------------------
# Function and its exact derivative 
f = function(x) x^3*sin(x/3)*log(sqrt(x))
dif_exact= function(x) (3*(x^2)*sin(x/3)*log(x)+(x^2)*sin(x/3))/2 + ((x^3)*cos(x/3)*log(x))/6

# Plot
x=seq(10,50,1)
dfplot<-data.frame(x=x,y=f(x),dif_exc=dif_exact(x)) %>% 
        gather(value,var,-x)
ggplot(dfplot,aes(x=x,y=var,color=value))+
  geom_line()+
  theme_classic()

# Error as a function of h
x0 = 1
err = numeric(16)
for (i in 1:16) {
  h = 10^-i
  err[i] = abs( (f(x0+h)-f(x0-h))/(2*h) - 0.5*sin(1/3) )
}
plot(log10(err), type="b", xlab="-log10(h)")

# numerical diff, richardson method
dif_exact(1)
err_s = abs(dif_exact(1) - grad(f,1,method="simple")) #fast method
err_r = abs(dif_exact(1) - grad(f,1,method="Richardson")) #accurate 
err_c = abs(dif_exact(1) - grad(f,1,method="complex")) #fast and accurate
# complex functions with the condition that x0 and f(x=0) are real

# Plot of numerical aproximation
x=seq(19,21,0.01)
dfplot2<-data.frame(x=x,y=f(x),dif_exc=dif_exact(x)) %>% 
         mutate(numdif_s=grad(f,x,method="simple")) %>% 
         gather(value,var,-x) 
ggplot(dfplot2,aes(x=x,y=var,color=value))+
  geom_line()+
  theme_classic()

# Gradiante
# Multivariate functions
f<- function(u){
  x = u[1]; y = u[2]; z = u[3]
  return(2*x + 3*y^2 - sin(z))
}
p_eval<-c(1,1,0)
grad(f,p_eval)

# Jacobian
F = function(x) c(x[1]^2 + 2*x[2]^2 - 3, cos(pi*x[1]/2) -5*x[2]^3)
jacobian(F, c(2,1))

# Hessian
hessian(f,c(1,1,0))

## Using the pracma package
f = function(x) x^3*sin(x/3)*log(sqrt(x))
x = 1:4
fderiv(f,x) # first derivative 
fderiv(f,x,n=2) # second derivative 

# numerical diff an error by Richardson's extrapolation
numderiv(f, x0=1, h=1/2) # for a escalar
numdiff(f,x=2:4) # for a vector

# this package allows us to draw the grid
v = seq(-2, 2, by=0.2)
X = meshgrid(v, v)$X
Y = meshgrid(v, v)$Y
Z = -(1/sqrt((X+1)^2 + Y^2) - 1/sqrt((X-1)^2 + Y^2))
par(mar=c(4,4,1.5,1.5),mex=.8,mgp=c(2,.5,0),tcl=0.3)
contour(v, v, t(Z), col="black",xlab="x",ylab="y")
grid(col="white")
grX = gradient(Z, v, v)$X
grY = gradient(Z, v, v)$Y
quiver(X, Y, grX, grY, scale = 0.2, col="black")

# Laplacian
# Differential operator that is the sum of the second derivatives
f = function(x) 2/x[1] - 1/x[2]^2
laplacian(f, c(1,1))

# 2.- Numerical Integration
#------------------------------
f = function(x) exp(-x)*cos(x)
curve(f,0,pi);grid()
q = integrate(f, 0, pi)
q

f = function(x) x^(1/3)/(1+x)
curve(f,0,1)
integrate(f,0,1)

# Integrate discretized functions
xs = seq(0, pi, length.out = 101)
ys = f(xs)
trapz(xs, ys) # high error
fsp = splinefun(xs, ys)
integrate(fsp, 0, pi)

# Gaussian Quadrature 
f = function(x) x^6
Lq = legendre.quadrature.rules(4)[[4]]
legendre.quadrature(f, Lq, lower = -1, upper = 1)

## exact value
intf = function(x) x^7/7
intf(1)-intf(-1)

## possible problems with oscillating functions
f = function(x) sin(1/x)
curve(f,0,1);grid()
integrate(f, 0, 1)
quadgk(f, 0, 1) #using gaussian cuadrature

## Multidimensional integration 
f = function(x) 1/(1 + x[1]^2 + x[2]^2)
q2 = adaptIntegrate(f, c(0, 0), c(1, 1))
q2
f = function(x, y) 1 / (1 + x^2 + y^2)
q3 = integral2(f, 0, 1, 0, 1)
q3

## 3.- Symbolic Manipulation in R
#---------------------------------
f = expression(sin(x)*exp(-a*x))
ffun = function(x,a) eval(f) 

# First Derivative 
g = D(f,"x")  #as symbolic
gfun = function(x,a) eval(g) # as function

# Second derivative
g2 = D(g,"x")
gfun2 = function(x,a) eval(g2) # as function

# Plot
curve(ffun(x,1),0,4, ylim = c(-1,1), ylab=c("f(x,1) and derivatives"))
curve(gfun(x,1), add=T, lty=2)
curve(gfun2(x,1), add=T, lty=3)
legend("topright", legend = c("f(x,1)", "df/dx", "d2f/dx2"),lty=1:3, bty="n")

# Other example
f2=expression(x^3-x^2+10)
ff2=function(x) eval(f2)
fd1=D(f2,"x")
gfun1=function(x) eval(fd1)
fd2=D(fd1,"x")
gfun2=function(x) eval(fd2)

curve(ff2(x),-3,3, ylim = c(0,20), ylab=c("f(x,1) and derivatives"))
curve(gfun1(x), add=T, lty=2)
curve(gfun2(x), add=T, lty=3)
legend("topright", legend = c("f(x)", "df/dx", "d2f/dx2"),lty=1:3, bty="n")

## 4.- Optimization
#--------------------
# 4.1.- One dimension
f = function(x) x*(20-2*x)*(16-2*x)
curve(f,0,8)
optimize(f, c(0,8), maximum=T)

f = function(x) x*sin(4*x)
curve(f,0,3)
optimize(f, c(0,3), maximum=T)
optimize(f, c(0,3)) # gives the first minimum encountered
optimize(f, c(1.5,3))

f.mins = findmins(f,0,3)
f.mins

f = function(x) abs(x^2-8)
curve(f,-4,4)
findmins(f,-4,4)

# 4.2 Multidimensional optimization
x1 = x2 = seq(.1,.9,.02)
z = outer(x1,x2,FUN=function(x1,x2) 1/x1 + 1/x2 +(1-x2)/(x2*(1-x1)) + 1/((1-x1)*(1-x2)))
persp(x1,x2,z,theta=45,phi=0)

f = function(x) {
  x1 = x[1]
  x2 = x[2]
  return(1/x1 + 1/x2 + (1-x2)/(x2*(1-x1)) + 1/((1-x1)*(1-x2)))
}
optim(c(.5,.5),f)

# a more accurate result is obtain if the gradient is given
fr = function(x) { # Rosenbrock Banana function
  x1 = x[1]
  x2 = x[2]
  100*(x2 - x1*x2)^2 + (1 - x1)^2
}
optim(c(-1.2,1), fr)

grr = function(x) { ## Gradient of fr
  x1 = x[1]
  x2 = x[2]
  c(-400*x1*(x2-x1*x1)-2*(1-x1),200*(x2 - x1*x1))
}
optim(c(-1.2,1), fr, grr, method = "BFGS")
optim(c(-1.2,1), fr, grr, method = "BFGS",hessian=TRUE)

# least-square fitting
set.seed(237)
x = seq(0, pi, length.out = 50)
y = sin(x) + 0.1*rnorm(50)
plot(x, y)

xp = seq(0, pi, length.out = 12)
F = function(p) {
  fsp = splinefun(xp, c(0, p, 0))
  sum((y - fsp(x))^2)
}
opt = optim(rep(0.5, 10), F, method="L-BFGS-B",
            lower = rep(0, 10), upper = rep(1, 10))

fsp = splinefun(xp, c(0, opt$par, 0))
yy = fsp(x)
lines(x, yy)

# Using nlm package
nlm(fr,c(-2,2))
f <- function(x, a) sum((x-a)^2)
nlm(f, c(10,10), a = c(3,5))

f <- function(x, a){
  res <- sum((x-a)^2)
  attr(res, "gradient") <- 2*(x-a)
  res
}
nlm(f, c(10,10), a = c(3,5))

# ucminf package
f1 = function(x) (2*x[1]-5)^2 + (x[2]-3)^2 + (5*x[3]-2)^2
ucminf(c(1,1,1), f1)

# the BB package
f2= function(x) -sin(x[1])*sin(x[2])*sin(x[1]+x[2])
spg(c(pi/2, pi/2), f2)

# 4.3.- Optimization with constrains
fr = function(x) {
  x1 = x[1]
  x2 = x[2]
  100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
grr <- function(x) { # gradient
  x1 = x[1]
  x2 = x[2]
  c(-400 * x1 * (x2 - x1 * x1) - 2 * (1 - x1), 200 * (x2 - x1 * x1))
}
constrOptim(c(-1.2,0.9), #intial gueses
            fr, #function
            grr, #gradient
            ui=rbind(c(-1,0),c(0,-1)), # constrains
            ci=c(-1,-1)) #length of the constrain 

# Alabama package
f = function(x) sin(x[1]*x[2]+x[3])
heq = function(x) -x[1]*x[2]^3+x[1]^2*x[3]^2-5 #equality constrain
hin = function(x) {
  h = rep(NA,2)
  h[1] = x[1]-x[2]
  h[2] = x[2] -x[3]
  h
  } #inequality constrain
p0 = c(3,2,1) #inital guess
constrOptim.nl(par=p0, fn = f, heq=heq, hin = hin)

# Rsolnp package
lower = rep(0,3)
upper = rep(5,3)
solnp(pars=p0, fun = f, eqfun=heq, ineqfun = hin, LB=lower,
      UB=upper, ineqLB = c(0,0), ineqUB = c(5,5))

# 4.4.- Linear Programming
obj = c(500, 400) # objective function 500*m1+400m2
mat = matrix(c(20, 20, # constraint matrix
               5, 30,
               15, 7),
      nrow=3, ncol=2, byrow=TRUE)
dir = rep("<=", 3) #direction of inequalities
rhs = c(100, 50, 60) # right hand side of inequalities
soln = lp("max", obj, mat, dir, rhs)
soln
soln$solution
mat %*% soln$solution

# Package linprog
solveLP(obj, rhs, mat, maximum = TRUE)

# 4.5.- Quadratic Programming
Dmat = matrix(c(1,-1,-1,2),2,2)
dvec = c(2,6)
Amat = matrix(c(-1,-1,1,-2,-2,-1),2,3)
bvec = c(-2,-2,-3)
solve.QP(Dmat, dvec, Amat, bvec)

# 4.6.- mixed-integer linear programs
obj = c(500, 450)
A = matrix(c( 6, 5,
              10, 20,
              1, 0), ncol = 2, byrow = TRUE)
b = c(60, 150, 8)
# Declare which variables have to be integer (here all of them)
int.vec = c(1, 2)
soln = lp("max", obj, A, rep("<=", 3), b, int.vec = int.vec)
soln
soln$solution

# 4.7.- Transportation problems
C = matrix(c(10, 70, 40,
             20, 30, 50), nrow = 2, byrow=TRUE)
row.dir = rep("<=", 2)
row.rhs = c(300, 200)
col.dir = rep(">=", 3)
col.rhs = c(150, 250, 100)
T = lp.transport(C, "min", row.dir, row.rhs, col.dir, col.rhs)
T
T$solution

# 4.8.- Assignment problems
M <- matrix(c(NA, 8, 6, 12, 1,
              15, 12, 7, NA, 10,
              10, NA, 5, 14, NA,
              12, NA, 12, 16, 15,
              18, 17, 14, NA, 13), 5, 5, byrow = TRUE)
M[is.na(M)] = 100 # replace NA with high value for constrain
A = lp.assign(M)
A
A$solution

# 4.9.- Subsetsum problems  GOOD
p = c(99.28, 5.79, 63.31, 89.36, 7.63, 30.77, 23.54, 84.24,
      93.29, 53.47, 88.19, 91.49, 34.46, 52.13, 43.09, 76.40,
      21.42, 63.64, 28.79, 73.03, 8.29, 92.06, 26.69, 89.07,
      10.03, 10.24, 40.29, 81.76, 49.01, 3.85) # possible values
P = as.integer(100*p)
obj = P
M = rbind(rep(1, 30), P)
dir = c("==", "<=")
rhs = c(4, 20010)
binary.vec = 1:40
L = lp("max", obj, M, dir, rhs, binary.vec = binary.vec)
L
inds = which(L$solution == 1)
inds; P[inds]/100; sum(P[inds])/100

## 4.- Ordinary Differential Equations
#--------------------------------------
# Compute parameter
G = 3600^2*6.673e-20
Msun = 1.989e30
GM = G*Msun
parms = GM

#Initialize variables
x0 = 149.6e6; vx0 = 0
y0 = 0; vy0 = 29.786*3600

# Set time sequence
tmin = 0; tmax = 8800; dt = 400
hrs = seq(tmin, tmax, dt)

# Set model
orbit = function(t,y,GM) {
 dy1 = y[2]
 dy2 = -GM*y[1]/(y[1]^2+y[3]^2)^1.5
 dy3 = y[4]
 dy4 = -GM*y[3]/(y[1]^2+y[3]^2)^1.5
 return(list(c(dy1,dy2,dy3,dy4)))
}

# Solve ode
out = ode(c(x0, vx0, y0, vy0), hrs, orbit, parms)
hrs = out[,1]; x = out[,2]; vx = out[,3]
y = out[,4]; vy = out[,5]
r = round(sqrt(x^2 + y^2)*1e-8,3)
v = round(sqrt(vx^2 + vy^2)/3600,3)
mat = cbind(hrs,x,y,r,v)
colnames(mat) = c("hrs", "x km", "y km",
                    "r/1e8 km", "v km/s")
mat

# Bessel differential equation
diffeqs = function(x,y,nu) {
  J=y[1]
  dJdx = y[2]
  with(as.list(parms), {
    dJ = dJdx
    ddJdx = -1/x^2*(x*dJdx + (x^2-nu^2)*J)
    res = c(dJ, ddJdx)
    list(res)
  })
}

# solve
 # Abscissa steps
xmin = 1e-15 # Don't start exactly at zero, to avoid infinity
xmax = 15
dx = 0.1
xx = seq(xmin, xmax, dx)

# Parameters
parms = c(nu = 1) # Bessel equation of order 1

# Initial values
y0 = c(J = 0, dJdx = 1)

# Solve with lsoda
out = lsoda(y0, xx, diffeqs, parms)

# Plot results and compare with built-in besselJ
xx = out[,1]; J = out[,2]; dJdx = out[,3]
plot(xx, J, type="l"); curve(besselJ(x,1),0,15,add=TRUE)
abline(0,0)

# Reaction mechanism
diffeqs = function(t,x,parms) {
  X=x[1]
  Y=x[2]
  with(as.list(parms), {
  dX = k1*A - k2*X - k3*X*Y^2
  dY = k2*X + k3*X*Y^2 - k4*Y
  list(c(dX, dY))
  })}

#times steps
tmin = 0; tmax = 200; dt = 0.01
times = seq(tmin, tmax, dt)

# Parameters: Rate constants and A concentration
parms = c(k1 = 0.01, k2 = 0.01, k3 = 1e6, k4 = 1, A = 0.05)

# Initial values
x0 = c(X = 0, Y = 0)

# Solve with adams method
out = ode(x0, times, diffeqs, parms, method = "adams")

# Plot results
time = out[,1]; X = out[,2]; Y = out[,3]
par(mfrow = c(1,3))
plot(time, X, type="l") # Time variation of X
plot(time, Y, type="l") # Time variation of Y
plot(X,Y, type="l") # Phase plot of X vs Y
par(mfrow = c(1,1))

# Population dynamics
diffeqs = function(t,x,parms) {
  prey = x[1]
  pred = x[2]
  with(as.list(parms), {
    dprey = -k1*pred*prey + k2*prey
    dpred = k3*pred*prey - k4*pred
    res = c(dprey, dpred)
    list(res)
    })
}

# Time steps
tmin = 0; tmax = 200; dt = 1
times = seq(tmin, tmax, dt)

# Parameters
parms = c(k1 = 0.01, k2 = 0.1, k3 = 0.001, k4 = 0.05)
# Initial values
x0 = c(prey = 50, pred = 15)

# Solve with lsoda
out = lsoda(x0, times, diffeqs, parms, rtol = 1e-6,
            atol = 1e-6)

# Plot results
par(mfrow = c(1,2))
time = out[,1]; prey = out[,2]; pred = out[,3]
plot(time, prey, type="l", ylim = c(0,100))
lines(time, pred, lty = 2)
legend("topleft", legend = c("prey","pred"), lty = 1:2,
       bty="n")
plot(prey, pred, type = "l")
par(mfrow = c(1,1))

# Exercise
# 1.- solve  dy/dt = ty+t initial condition y0=0.5
ode1<-function(t,y,parm){
  y=Y[1]
  dydt = t*y+t
  list(dydt)
}
yinit = 0.5
times = seq(0,5,0.1)
out = as.data.frame(ode(y = yinit, 
                        func = ode1, 
                        times, 
                        parms = c()))
names(out) = c("t","y")
out
plot(out$t,out$y,type="l")

## 5.- Fitting models to data
#-----------------------------
data<-read.table(file="data.txt")
data<-data[-1,] %>% select("V2","V3") %>% as.data.frame()
names(data)<-c("y","x")
x=data$x %>% as.numeric();  y=data$y %>% as.numeric()
nlsfit2 = nls(y~b1*(1-exp(-b2*x)), start=list(b1=250,b2=5e-4))
summary(nlsfit2)
nlsLM_lan3 = nlsLM(y~b1*exp(-b2*x)+b3*exp(-b4*x)+b5*exp(-+ b6*x),
                   start=list(b1=1.2,b2=0.3,b3=5.6,b4=5.5,b5=6.5,b6=7.6),
                  control = nls.lm.control(maxiter=100))
summary(nlsLM_lan3)

data <- read_xlsx("iphone_sales.xlsx") %>% as.data.frame()
year_sales <- data %>% 
  group_by(Year) %>% 
  mutate(sale=sum(Sales)) %>%
  select("Year","sale") %>% 
  distinct(Year, sale,.keep_all = TRUE) %>% 
  as.data.frame()
Year<-year_sales$Year
sales <-year_sales$sale
t<-1:nrow(year_sales)
bass_model<-formula(sales ~ M*(((P+Q)^2/P)*exp(-(P+Q)*t))/((1+(Q/P)*exp(-(P+Q)*t))^2))

## Gauss-Newton
gn_lm <-nls(formula = bass_model,
            start = c(M=1000,P=0.01,Q=0.5),
            trace = TRUE)
gn_lm

### Lenvember Marquart
LM_algorithm<-nlsLM(bass_model,
                    start=list(M=1000,P=0.01,Q=0.5),
                    algorithm = "LM",
                    trace=TRUE) 
summary(LM_algorithm)  

