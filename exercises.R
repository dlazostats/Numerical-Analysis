# Exercises in Elements of Numerical Analysis  
#============================================
library(deSolve)
library(dplyr)
library(ggplot2)
options(scipen=999)

# Working directory
script_name <- 'exercises.R'
ruta <- gsub(rstudioapi::getActiveDocumentContext()$path,pattern = script_name,replacement = '')
setwd(ruta)
source("euler.R")

# Exercise 1
ode1<-function(t,y,params){
  dydt = 1-t+4*y
  list(dydt)
}
sol_ex1<-function(t){
  y=(19/16)*exp(4*t)+(1/4)*t-(3/16)
  return(y)
}
yini<-1
times <- seq(from = 0, to = 1, by = 0.1)
s11 <- ode(func = ode1,y = yini, times = times,
          parms = NULL,method="euler")
s12 <- ode(func = ode1,y = yini, times = times,
          parms = NULL,method="rk4")
s13 <- ode(func = ode1,y = yini, times = times,
           parms = NULL,method="adams")
sol_exdf1<-data.frame(time=times,y=sol_ex1(times),method="Exact")
dfp<-rbind(s11,s12) %>% 
     as.data.frame() %>% 
     mutate(method=c(rep("Euler",nrow(s11)),rep("Rk4",nrow(s11))))
names(dfp)<-c("time","y","method")
ggplot(rbind(dfp,sol_exdf1),aes(x=time,y=y,color=method))+
  geom_line()+
  geom_point()+
  theme_classic()

## error
df_er<-cbind(s11,s12,s13,sol_exdf1) 
df_er<-df_er[,c(1,2,4,6,8)]
names(df_er)<-c("time","y_eule","y_rk4","y_adams","y_exact")
df_er$err_eul<-abs(df_er$y_eule-df_er$y_exact)
df_er$err_rk<-abs(df_er$y_rk4-df_er$y_exact)
df_er$err_adm<-abs(df_er$y_adams-df_er$y_exact)

