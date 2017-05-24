


plot(function(x) dgamma(x,shape=3.2,rate=3.2),xlim=c(0,3))

fgammamod = function(x,L, n) {
    
    (x**(n-1))*exp(-L*x)    *  (L**n)/gamma(n) 
}

#plot(function(x) f(x,1,3.2),xlim=c(0,3))


fstar = function(x,q,l) {
    
    
    sum1 = 0
    for(k in 1:5) 
        sum1 = sum1 + fgammamod(x, 2*q*(l+1) , k*(l+1))/ (2**k)
    
    sum1
    
}


plot(function(x) fstar(x,1,4.5),xlim=c(0,3), xlab="Morgans", ylab="Density")
title("Crossover distance")
abline(h=1,col="gray")





dens = sapply(seq(0,4,by=0.01), function(x) fstar(x,1,4.5) )

cdf = cumsum(dens)
cdf = cdf/max(cdf)

cdf[1:10]
crossoverCDFvector = cdf


generateRecombinationDistances = function ( n ) {
    
    t = findInterval( runif(n), crossoverCDFvector )
    t = t * 4 / length(crossoverCDFvector)
    100 * t
}

generateRecombinationDistances(20)

t=findInterval(runif(100000),cdf)
hist(t,n=100)    

mean(t)






t=rexp(100000,3)+rexp(100000,3)+rexp(100000,3)#+rexp(100000,3)++rexp(100000,3)
hist(t,n=200)
mean(t)

n=1e6
p=2*(1+4)
t = rchisq(n,p) + rchisq(n,p)+ rchisq(n,p)+ rchisq(n,p)+ rchisq(n,p)
t = t / mean(t)

hist(t,n=200,xlim=c(0,3))
