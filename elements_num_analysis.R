# Elements of Numerical Analysis  
#===============================
# 1.- Non Linear Algebraic Equations
#===================================

# 1.1.- Bisection method
#-----------------------
bisection<-function(f,a,b,tol=1e-3,m=100){
  require(dplyr)
  require(ggplot2)
  aini<-a;bini<-b
  iter<-0
  fa<-f(a)
  fb<-f(b)
  l<-list()
  while(abs(b-a)>tol && (iter < m)){
    iter <- iter + 1
    if (iter > m){
       warning (" iterations maximum exceeded ")
       break
    } 
    xmid <- (a + b)/2
    ymid <- f(xmid)
    if(fa*ymid > 0){
      a <- xmid
      fa <- ymid
    }else{
      b <- xmid
      fb <- ymid
    }
    error<-xmid-b
    l[[iter]]<-c(iter,a,b,xmid,ymid,error)
  }
  dfr<-do.call("rbind",l) %>% as.data.frame()
  names(dfr)<-c("n","a","b","xmind","ymind","tol")
  root<-(a+b)/2
  lr<-list()
  lr$root<-root
  lr$table<-dfr
  lr$iter<-nrow(dfr)
  #lr
  xp<-seq(aini,bini,(bini-aini)/100)
  dfp<-data.frame(x=xp,y=f(xp)) 
  p1<-ggplot(dfp,aes(x=x,y=y))+
    geom_line()+
    geom_hline(yintercept=0, linetype="dashed", color = "red")+
    theme_classic()
  lr$plot<-p1
  lr
}

# Example 1
f<-function(x){x^2-1}
root1<-bisection(f,0.5,1.25)
root1$root
root1$plot

# Example 2
f<-function(x){exp(x)-sin(x)}
root2<-bisection(f,-8,-4)
root2$root
root2$plot
root2$table

# 1.2.- Fixed Point
#---------------------
fixedpoint <- function(fun, x0, tol=1e-07, niter=100){
  require(dplyr)
  xold <- x0
  xnew <- fun(xold)
  iter <- 1
  l<-list();lr<-list()
  while ((abs(xnew-xold) > tol) && (iter < niter)) {
    xold <- xnew;
    xnew <- fun(xold);
    iter <- iter + 1
    if( iter > niter)
      stop ("No solution found ")
    error<-abs(xnew-xold)
    l[[iter]]<-c(iter,xnew,error)
  }
  dfr<-do.call("rbind",l) %>% as.data.frame()
  names(dfr)<-c("n","x","tol")
  lr$root<-xnew
  lr$Iter<-nrow(dfr)
  lr$tabla<-dfr
  lr
}

#Example 1
f<-function(x){x - log(x) + exp(-x)}
fixedpoint(f,x=2)

# 1.3.- Newton-Raphson
#---------------------
newton_raph <- function(fun, x, tol = 1e-3, m = 100) {
  iter <- 0
  oldx <- x
  x <- oldx + 10*tol
  g = D(fun,"x")  
  fp = function(x) eval(g)
  l<-list()
  f <- function(x){eval(fun[[1]])}
  while (abs(x - oldx) > tol && (iter < m)) {
    iter <- iter + 1
    if( iter > m)
      stop ("No solution found ")
    oldx <- x
    x <- x - f(x)/fp(x)
    error <- abs(oldx-x)
    l[[iter]]<-c(iter,x,error)
  }
  dfr<-do.call("rbind",l) %>% as.data.frame()
  names(dfr)<-c("n","x","tol")
  root<-x
  lr<-list()
  lr$root<-root
  lr$Iter<-nrow(dfr)
  lr$tabla<-dfr
  lr$deriv<-g
  lr
}

#Example 1
f<-expression(x^2 - 2*x + 1)
newton_raph(f,x=1.25)

# 1.4.- Newton-Raphson-Modificated
#----------------------------------
newton_raph_mod <- function(fun, x, tol = 1e-3, m = 100) {
  f <- function(x){eval(fun[[1]])}
  iter <- 0
  oldx <- x
  x <- oldx + 10*tol
  g1 = D(fun,"x")  
  fp1 = function(x) eval(g1)
  g2 = D(g1,"x")  
  fp2 = function(x) eval(g2)
  l<-list()
  while (abs(x - oldx) > tol && (iter < m)) {
    iter <- iter + 1
    if( iter > m)
      stop ("No solution found ")
    oldx <- x
    x <- x - (fp1(x)*f(x))/(fp1(x)^2-fp2(x)*f(x))
    error <- abs(oldx-x)
    l[[iter]]<-c(iter,x,error)
  }
  dfr<-do.call("rbind",l) %>% as.data.frame()
  names(dfr)<-c("n","x","tol")
  root<-x
  lr<-list()
  lr$root<-root
  lr$Iter<-nrow(dfr)
  lr$tabla<-dfr
  lr
}

#Example 1
f<-expression(10*x^3 - 8.3*x^2 +2.295*x-0.21141)
newton_raph_mod(f,x=0.265)

# 1.5.- Secant Method
#---------------------
secant <- function (f, x, tol = 1e-3, m = 100) {
  i <- 0
  oldx <- x
  oldfx <- f(x)
  x <- oldx + 10*tol
  l<-list();lr<-list()
  while (abs(x - oldx ) > tol  && (i < m)) {
    i <- i + 1
    if (i > m)
      stop ("No solution found ")
    fx <- f(x)
    newx <- x - fx*((x - oldx )/(fx - oldfx ))
    oldx <- x
    oldfx <- fx
    x <- newx
    error<-abs(x-oldx)
    l[[i]]<-c(i,x,error)
  }
  dfr<-do.call("rbind",l) %>% as.data.frame()
  names(dfr)<-c("n","x","tol")
  lr$root<-x
  lr$Iter<-nrow(dfr)
  lr$tabla<-dfr
  lr
}
f<-function(x){x^2 - 1}
secant(f,x=50)

# 2 .- Optimization
#===================
# 2.1.- Golden Search
#---------------------
goldsectmax <- function (f, a, b, tol = 1e-3, m = 100) {
  iter <- 0
  phi <- ( sqrt (5) - 1) / 2
  a.star <- b - phi*abs(b - a)
  b.star <- a + phi*abs(b - a)
  l<-list();lr<-list()
  while (abs(b - a) > tol && (iter<m)) {
    iter <- iter + 1
    if ( iter > m) {
      warning (" iterations maximum exceeded ")
      break
    }
    if(f(a.star ) > f(b.star )) {
      b <- b.star
      b.star <- a.star
      a.star <- b - phi*abs(b - a)
    } else {
      a <- a.star
      a.star <- b.star
      b.star <- a + phi*abs(b - a)
    }
    xopt<-(a+b)/2
    error<-abs(b-a)
    l[[iter]]<-c(iter,a,b,xopt,error)
  }
  dfr<-do.call("rbind",l) %>% as.data.frame()
  names(dfr)<-c("n","a","b","x","tol")
  lr$root<-(a + b)/2
  lr$Iter<-nrow(dfr)
  lr$tabla<-dfr
  lr
}

#Example 1
f <- function(x) { -x^2 + 4*x + 1 }
goldsectmax(f, 1, 3)
xx<-seq(1,3,0.1)
y<-f(xx)
plot(xx,y,type="l")

#Example 2
f <- function(x) { sin(x) - x^2 }
goldsectmax(f, 0, 1)

#Example 2
fp <- function(x){(1/4)*x^4 + x^3 - x - 1 }
xx<-seq(-1,1,0.1)
y<-fp(xx)
plot(xx,y,type="l")

# max
goldsectmax(fp, -1, 1)

# min
fpm <- function(x){-((1/4)*x^4 + x^3 - x - 1)}
goldsectmax(fpm, -1, 1)

# 2.2.- Gradient Descent
#-----------------------
graddsc <- function (fun , x, h = 1e-3, tol = 1e-4, m = 1e3){
  g = D(fun,"x")  
  fp = function(x) eval(g)
  iter <- 0
  oldx <- x
  x = x - h*fp(x)
  l<-list();lr<-list()
  while(abs(x - oldx) > tol && (iter<m)) {
    iter <- iter + 1
    if( iter > m)
      stop ("No solution found ")
    oldx <- x
    x = x - h*fp(x)
    xopt<-x
    error<-abs(x-oldx)
    l[[iter]]<-c(iter,xopt,error)
  }
  dfr<-do.call("rbind",l) %>% as.data.frame()
  names(dfr)<-c("n","x","tol")
  lr$root<-x
  lr$Iter<-nrow(dfr)
  lr$tabla<-dfr
  lr$deriv<-g
  lr
}
fp <- expression((1/4)*x^4 + x^3 - x - 1)
graddsc(fp, x = -2, h = 1/100)

# 2.3.- Multidimensional Gradient Descent
#----------------------------------------
mul_graddsc<- function (fp , x, h = 1e2 , tol = 1e-4, m = 1e3) {
  require(splus2R)
  iter <- 0
  oldx <- x
  x = x - h*fp(x)
  l<-list();lr<-list()
  while ( vecnorm(x - oldx) > tol && (iter<m)){
    iter <- iter + 1
    if( iter > m)
      return("No solution found ")
    oldx <- x
    x = x - h*fp(x)
    xopt<-x
    error<-abs(x-oldx)
    l[[iter]]<-c(iter,xopt,error)
  }
  dfr<-do.call("rbind",l) %>% as.data.frame()
  names(dfr)<-c("n","x1","x2","tol")
  lr$root<-x
  lr$Iter<-nrow(dfr)
  lr$tabla<-dfr
  lr
}
fp<-function(x){
  x1 = 2*x[1] - 2
  x2 = 8*x[2] - 8
  return(c(x1,x2))
}
mul_graddsc(fp,c(0,0),0.05,m=100)

# 2.4.- Hill Climbing 
#--------------------
hillclimbing <- function (f, x, h = 1, m = 1e3) {
  #set.seed(123)
  n <- length(x)
  xcurr <- x
  ycurr <- f(x)
  l<-list();lr<-list()
  for(i in 1:m) {
    xnext <- xcurr
    i <- ceiling( runif(1, 0, n))
    xnext [i] <- rnorm (1, xcurr[i], h)
    ynext <- f( xnext )
    if( ynext < ycurr ) {
      xcurr <- xnext
      ycurr <- ynext
    }
    xopt<-xcurr
    l[[i]]<-c(i,xopt)
  }
  dfr<-do.call("rbind",l) %>% as.data.frame()
  names(dfr)<-c("n",paste0("var",seq(1,ncol(dfr)-1)))
  lr$root<-xcurr
  lr$Iter<-nrow(dfr)
  lr$tabla<-dfr
  lr
  #return(xcurr)
}

#Example 1
f <- function(x){(x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2}
hillclimbing(f,c(0,0))

#Example 2
f <- function(x){
  (x[1]*x[2]*x[3])
}
hillclimbing(f,c(0.5,0.5,0.5))

# using packages
f1 <- function(x) 2*(x[1]-1)^2 + 5*(x[2]-3)^2 +10
r1 <- optim(c(1, 1), f1)
r1

f2 <- function(x){(x[1]^2 + x[2] - 11)^2 + (x[1] + x[2]^2 - 7)^2}
r2 <- optim(c(1, 1), f2)
r2

f <- function(x)(print(x) - 1/3)^2
xmin <- optimize(f,
                 interval = c(0, 1),
                 tol = 0.0001)
xmin

# 3 .- Differential Equations
#============================
# 3.1 Euler Method
#-----------------
euler <- function (f, x0 , y0 , h, n) {
  x <- x0
  y <- y0
  for(i in 1:n) {
    y0 <- y0 + h*f(x0 , y0)
    x0 <- x0 + h
    x <- c(x, x0)
    y <- c(y, y0)
  }
  return ( data.frame (x = x, y = y))
}

#Example 1
f <- function(x, y) { y/(2*x + 1)}
euler(f, 0, 1, h=0.5,n=5)

#Example 2
f<-function(x,y){ -x*y^2}
euler(f, 0, 2, h=0.1,n=10)

# 3.2 Midpoint
#--------------
midptivp <- function(f, x0 , y0 , h, n) {
  x <- x0
  y <- y0
  for(i in 1:n) {
    s1 <- h*f(x0 , y0)
    s2 <- h*f(x0 + h / 2, y0+s1/2)
    y0 <- y0 + s2
    x0 <- x0 + h
    x <- c(x, x0)
    y <- c(y, y0)
  }
  return ( data.frame (x = x, y = y))
}
f <- function(x, y) { y/(2*x+1)}
midptivp(f, 0, 1, 0.5, 5)

# 3.3 Runge-Kutta
#----------------
rungekutta4 <- function (f, x0 , y0 , h, n) {
  x <- x0
  y <- y0
  for(i in 1:n) {
    s1 <- h*f(x0 , y0)
    s2 <- h*f(x0 + h/2, y0 + s1/2)
    s3 <- h*f(x0 + h/2, y0 + s2/2)
    s4 <- h*f(x0 + h, y0 + s3)
    y0 <- y0 + s1/6 + s2/3 + s3/3 + s4/6
    x0 <- x0 + h
    x <- c(x, x0)
    y <- c(y, y0)
  }
  return (data.frame(x = x, y = y))
}
f <- function(x, y) { y/(2*x+1)}
rungekutta4(f, 0, 1, 0.5, 5)

library(pracma)
rk4(f, 0, 1, 1, 5)

# 3.4 Adams Bashforth method
#----------------------------
adamsbashforth <- function (f, x0 , y0 , h, n) {
  ## Quick Euler the value of x1 , y1
  y1 <- y0 + h*f(x0 , y0)
  x1 <- x0 + h
  x <- c(x0 , x1)
  y <- c(y0 , y1)
  n <- n - 1
  for(i in 1:n) {
    yn <- y1 + 1.5*h*f(x1 , y1)-.5*h*f(x0 , y0)
    xn <- x1 + h
    y0 <- y1
    x0 <- x1
    y1 <- yn
    x1 <- xn
    y <- c(y, y1)
    x <- c(x, x1)
  }
  return ( data.frame (x = x, y = y))
}
f <- function(x, y) { y/(2*x+1)}
adamsbashforth(f, 0, 1, 0.5, 5)
