rm(list=ls())

#parameters
Niterations = 1000
Msamples.min = 5
Msamples.max = 100
sample.freq = 5
truechange = 1
value_nullhypothesis = 1

# define log likelihood ratio function
calc.llr <- function(r,len.x,len.y){
  len.x * log(len.x/len.y + 1/r) + len.y*log(len.y/len.x + r) + len.x*log(len.y/(len.x+len.y)) + len.y*log(len.x/(len.x+len.y))
}

#initialize
i = 0; j=0;
null_avg.exp = array(0, c(1,Niterations))
alt_avg.exp = array(0, c(1,Niterations))
ratio.means.exp = array(0, c(1,Niterations))
ratio.means.hi = array(0, c(1,Niterations))
likelihoodRatio = array(0, c(1,Niterations))
llr.obs.exp = array(0, c(1,Niterations))
p.value.exp = array(0, c(1,Niterations))
p.value.lo = array(0, c(1,Niterations))
p.value.hi = array(0, c(1,Niterations))
sig.exp = array(0, c(1,Niterations))
Pi.exp = array(0, c(1, 1 + (Msamples.max-Msamples.min) / sample.freq))
null.rate.pois = array(0, c(1,Niterations))
alt.rate.pois = array(0, c(1,Niterations))
null.ll.pois = array(0, c(1,Niterations))
alt.ll.pois = array(0, c(1,Niterations))
p.value.pois = array(0, c(1,Niterations))
sig.pois = array(0, c(1,Niterations))
Pi.pois = array(0, c(1, 1 + (Msamples.max-Msamples.min) / sample.freq))

for (Msamples.null in seq(Msamples.min, Msamples.max, sample.freq)){
    j = j+1
    Msamples.alt = Msamples.null * 1
    nullsamples.exp = array(0, c(Niterations, Msamples.null))
    alternativesamples.exp = array(0, c(Niterations, Msamples.alt))
    nullsamples.pois = nullsamples.exp
    alternativesamples.pois = alternativesamples.exp
    
    for (i in 1:Niterations){
      
      ######### Exponential
      nullsamples.exp[i,] = rexp(Msamples.null, 1/value_nullhypothesis)
      alternativesamples.exp[i,] = rexp(Msamples.alt, 1/(value_nullhypothesis+truechange))
      
      null_avg.exp[i] = sum(nullsamples.exp[i,]) / Msamples.null
      alt_avg.exp[i] = sum(alternativesamples.exp[i,]) / Msamples.alt
      
      # observed ratio of sample means
      ratio.means.exp[i] <- mean(nullsamples.exp[i,])/mean(alternativesamples.exp[i,])
      
      # observed log likelihood ratio
      llr.obs.exp[i] <- calc.llr(ratio.means.exp[i], Msamples.null, Msamples.alt) 
      
      # p-value in lower tail
      p.value.lo[i] <- pf(ratio.means.exp[i] ,2*Msamples.null,2*Msamples.alt) 
      
      # find the other ratio of sample means giving an LLR equal to that observed
      ratio.means.hi[i] <- uniroot(function(x) calc.llr(x,Msamples.null,Msamples.alt)-llr.obs.exp[i], lower=1, upper=100, tol=1e-6)$root
      
      #p.value in upper tail
      p.value.hi[i] <- 1-pf(ratio.means.hi[i], 2*Msamples.null, 2*Msamples.alt) 
      
      # overall p.value
      p.value.exp[i] <- p.value.lo[i] + p.value.hi[i]
      
      # sanity check
      p.value[i] <- 1-pchisq(2*llr.obs.exp[i], 1)
      
      # stat test two tailed
      sig.exp[i] <- p.value.exp[i]<0.05
      
      
      
      ######### Poisson
      nullsamples.pois[i,] = rpois(Msamples.null, value_nullhypothesis)
      alternativesamples.pois[i,] = rexp(Msamples.alt, value_nullhypothesis+truechange)
      
      null.rate.pois[i] = sum(nullsamples.pois[i,]) / Msamples.null
      alt.rate.pois[i] = sum(alternativesamples.pois[i,]) / Msamples.alt
      
      null.ll.pois[i] = sum(dpois(nullsamples.pois[i,], null.rate.pois[i], log=T))
      alt.ll.pois[i] = sum(dpois(nullsamples.pois[i,], alt.rate.pois[i], log=T))
      
      p.value.pois[i] = 1 - pchisq(2* (null.ll.pois[i] - alt.ll.pois[i]), df=1)
      sig.pois[i] = p.value.pois[i]<0.05
      
    }
    
    # Fraction of statistically significant results (power)
    Pi.exp[j] <- sum(sig.exp) / Niterations
    Pi.pois[j] <- sum(sig.pois) / Niterations
}

### PLOTTING POWER FUNCTIONS
png(file="powercalculations_saved_test.png",width=1800,height=900, res=300)

plot(seq(Msamples.min, Msamples.max, sample.freq), Pi.pois, col="red", xlim=c(0,100), ylim=c(0,1), type="b", xlab="",ylab="", axes=F)
par(new=T)
plot(seq(Msamples.min, Msamples.max, sample.freq), Pi.exp, col="blue", xlim=c(0,100), ylim=c(0,1), type="b",xlab="Sample size, X", ylab="Power", cex.lab=1.2)

legend(60, 0.6, c("Poisson", "Exponential"), lty=1, pch=1, lwd=1, col=c("red","blue"),bty="n")

c.exp = xy.coords(c(Pi.exp[6], Pi.exp[7]),c(6*sample.freq,7*sample.freq))
c.pois = xy.coords(c(Pi.pois[5], Pi.pois[6]),c(5*sample.freq,6*sample.freq))
xx.exp = approx(c.exp, xout=0.8, method="linear")
xx.pois = approx(c.pois, xout=0.8, method="linear")

par(new=T)
lines(c(xx.exp$y, xx.exp$y), c(0, 0.8), lty=2, col="blue")
par(new=T)
lines(c(xx.pois$y, xx.pois$y), c(0, 0.8), lty=2, col="red")

par(new=T)
lines(c(0, xx.exp$y), c(0.8, 0.8), lty=2, col="blue")
par(new=T)
lines(c(0, xx.pois$y), c(0.8, 0.8), lty=2, col="red")

dev.off()
