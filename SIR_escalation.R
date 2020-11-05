rm(list=ls())
setwd('~/Desktop/HW-re-escalation/code/')
require(deSolve)

#### Miller, Metcalfe, Grenfell model in Nature Medicine
bC <- 0.5
r<-0.6
bA <- 0.5
deltaE <- 1/3
deltaP <- 1/2.1
deltaC <- 1/2.9
deltaA <- 1/5
N<-1
tot.popn <-100000
# the times for the numerical integration
times <- seq(0, 150, by = .1)
# doubling time 7 days
lambda1 <- log(2)/7
# halving time 10 days.
lambda2 <- log(1/2)/10

# Given a user supplied lambda, calculate the beta
find.beta = function(beta,lambda){
  # Calculate R0 using the next-generation approach of van den Driesshe and
  # Watmough 2002.
  F <- matrix(rep(0,16), nrow =4, ncol = 4)
  V <- F
  F[1,] <- c(0, beta, beta*bC, beta*bA)
  V[1,1] <- deltaE
  V[2,1] <- -r*deltaE
  V[2,2] <- deltaP
  V[3,2] <- -deltaP
  V[3,3] <- deltaC
  V[4,1] <- -(1-r)*deltaE
  V[4,4] <- deltaA
  J <- F - V
  lambda.calc <- max(Re(eigen(J)$values))
  #inv.V = inv(V)
  #K = F%*%inv.V
  #R0.calc = max(abs(eigen(K)$values))
  #R0 - R0.calc
  lambda - lambda.calc
}

# Root finding to get the beta corresponding to the target lambda
# Some of the code in this section is copy & paste from the Miller et al. code
beta2<-uniroot(find.beta,c(0,10),lambda=lambda2)$root
beta1<-uniroot(find.beta,c(0,10),lambda=lambda1)$root

# beta, is assumed to have a breakpoint around t = SDstart which is the
# time that restrictions are enacted.
beta.fun = function(t){
  if(t>=t.hat){
    betaval <- beta2
  } else{
    betaval <- beta1}
  return(betaval)
}

Miller.CT = function(t,y,parms){
  S <- y[1]
  E <- y[2]
  IP <- y[3]
  IC <- y[4]
  IA <- y[5]
  cumI <- y[6]

  dS = -S*beta.fun(t)*(IP+bC*IC+bA*IA)/N
  dE = S*beta.fun(t)*(IP+bC*IC+bA*IA)/N-deltaE*E
  dIP = r*deltaE*E-deltaP*IP
  dIC = deltaP*IP-deltaC*IC
  dIA = (1-r)*deltaE*E - deltaA*IA
  cumI = S*beta.fun(t)*(IP+bC*IC+bA*IA)/N
  
  return(list(c(dS,dE,dIP,dIC,dIA,cumI)))
}

# Would like the initial condition to be on an eigenvector
t.hat <- 1000
eigen.vec = function(){
  I0=1/tot.popn
  # Initial conditions
  yini  = c(S = 1-I0, E = 0, IP = 0, IC = I0*r, IA = I0*(1-r), cumI=0)
  # performing the numerical integration
  out <- ode(y = yini, parms = NULL, times = times, func = Miller.CT)
  out<-data.frame(out)
  E = out$E[length(times)]
  IA = out$IA[length(times)]
  IP = out$IP[length(times)]
  IC = out$IC[length(times)]
  tot = E+IA+IP+IC
  E0 = E/tot
  IA0 = IA/tot
  IP0 = IP/tot
  IC0 = IC/tot
  return(data.frame(E0=E0, IP0 = IP0, IC0 = IC0, IA0=IA0))
}

# SIR model
gamma = 1/8
I0 <- 10/tot.popn

## Example figure
t.hat <- 30
tvec1 = seq(0,t.hat,.1)
tvec2 = seq(t.hat,100,.1)

png(file = "Results.png", width = 600)
par(mfrow = c(2,2), mar = c(5, 4, 2, 2))
plot(tvec1, I0*exp(lambda1*tvec1)*tot.popn, lty = 2, typ="l",xlim = c(0,100), ylab ="Number of infected individuals", xlab = "day", bty = "L")
lines(tvec2, I0*exp(lambda1*t.hat)*exp(lambda2*(tvec2-t.hat))*tot.popn, lty = 2)

ini = eigen.vec()
yini  = c(S = 1-I0, E = I0*ini$E0, IP = I0*ini$IP0, IC = I0*ini$IC0, IA = I0*ini$IA0, cumI = 0)
out <- ode(y = yini, parms = NULL, times = times, func = Miller.CT)
out <- data.frame(out)
# Number infectious
lines(out$time, (out$IA+out$IC+out$IP+out$E)*tot.popn, col = "black")

# Loop across day that restrictions enacted
t.hat.vec <- seq(1,45,1)

# Time at peak
t.hat = t.hat.vec
# Active infecteds at peak
I.max.SIR = I0*exp(lambda1*t.hat)
# Final size
final.size.SIR = I.max.SIR*(beta1/lambda1 + beta2/lambda2) - beta1*I0/lambda1
# Restrictions duration
tau.SIR = log(I0/I.max.SIR)/lambda2

# Preallocate output vector
out.data = matrix(0,nrow = length(t.hat.vec), ncol = length(times))
final.size = NULL
tau = NULL
peak.cases = NULL
Ic.data = out.data
tau.ic = NULL

# The is a loop across different days that escalation could start
for(i in seq(1,length(t.hat.vec))){
  t.hat <- t.hat.vec[i]
  # performing the numerical integration
  out <- ode(y = yini, parms = NULL, times = times, func = Miller.CT)
  out <- data.frame(out)
  out.data[i,]<-out$IC+out$E+out$IA+out$IP
  peak.cases[i] <-max(out$IC+out$E+out$IA+out$IP)*tot.popn
  final.size[i]<- out$cumI[length(times)]*tot.popn
  tau[i] <- out$time[max(which((out$IA+out$IC+out$IP+out$E)*tot.popn >= I0*tot.popn))] - t.hat
  Ic.data[i,] <- out$IC
}

# includes asymptomatic
plot(t.hat.vec, final.size, typ = "l", ylab = "total cases", xlab = "day restrictions enacted", bty = "L")

lines(t.hat.vec, final.size.SIR*tot.popn, lty = 2)

plot(t.hat.vec, tau, typ = "l", xlab = "day restrictions enacted", ylab = "duration of restrictions (days)", bty = "L")
lines(t.hat.vec, tau.SIR, lty=2)

plot(t.hat.vec, I.max.SIR*tot.popn, typ = "l", lty = 2, xlab = "day restrictions enacted", ylab = "peak number of infected individuals", bty = "L")
lines(t.hat.vec, peak.cases)
legend("topleft",
       legend = c("SIR", "COVID-19"), lty = c(2,1), bty = "n")
dev.off()

burn.in <- 30
scen1 <- burn.in + 7
scen2 <- burn.in + 14

i1 = times[max(which(Ic.data[scen1,]*tot.popn>=30))]-burn.in
i2 = times[max(which(Ic.data[scen2,]*tot.popn>=30))]-burn.in
t1 = max(which(times<=scen1))
t2 = max(which(times<=scen2))
 
png(file="two-scenarios.png",width=1100,height=550)
par(mfrow = c(1,2), mar = c(5.1, 5.1, 4.1, 2.1))
plot(times - burn.in, Ic.data[scen1,]*tot.popn, xlab="days", typ = "l",bty="L", ylab="number of clinical infections", xlim = c(0,50), main = "Restrictions enacted on day 7",ylim = c(0,150), lwd=3,cex.lab=2, cex.main = 2, cex.axis=2)
rect(scen1-burn.in, 0, i1, 30, col="orange", lwd=0)
lines(c(scen1-burn.in,scen1-burn.in),c(0,Ic.data[scen1,t1])*tot.popn, col = "orange", lwd = 3)
text(16,15,labels="restrictions", col = "white",cex=1.6)
text(15,90,labels="peak occurs after", col = "black", cex=1.5)
text(15,85,labels="restrictions enacted", col = "black", cex=1.5)
text(9,45,labels="day 7", col = "orange", cex=1.7, srt=90)

plot(times-burn.in, Ic.data[scen2,]*tot.popn, xlab="days", typ = "l",bty="L", ylab="number of clinical infections", xlim = c(0,50), main = "Restrictions enacted on day 14",ylim = c(0,150), lwd=3,cex.lab=2, cex.main = 2, cex.axis=2)
plot(times - burn.in, i2, xlab="days", typ = "l",bty="L", ylab="number of clinical infections", xlim = c(0,50), main = "Restrictions enacted on day 7",ylim = c(0,150), lwd=3,cex.lab=2, cex.main = 2, cex.axis=2)
rect(scen2-burn.in, 0, i2, 30, col="orange", lwd=0)
lines(c(scen2-burn.in,scen2-burn.in),c(0,Ic.data[scen2,t2])*tot.popn, col = "orange", lwd = 3)
text(27,20,labels="duration of restrictions", col = "white",cex=1.6)
text(29,10,labels="is longer", col = "white",cex=1.6)
text(16,50,labels="day 14", col = "orange", cex=1.7, srt=90)
dev.off()