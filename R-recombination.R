#### Standard Functions ####


a = read.table("geneticmap/genetic_map_GRCh37_chr4.txt",as=T,h=T)
colnames(a) = c("chr","bp","rate.cm","cm")
a = a[ order(a$cm),]

geneticMap = a




generateRecombinationDistances = function ( n ) {
    
    ( rexp(n, 0.01))
}





fgammamod = function(x,L, n) {
    
    (x**(n-1))*exp(-L*x)    *  (L**n)/gamma(n) 
}


fstar = function(x,q,l) {
    
    
    sum1 = 0
    for(k in 1:5) 
        sum1 = sum1 + fgammamod(x, 2*q*(l+1) , k*(l+1))/ (2**k)
    
    sum1
    
}





plot(function(x) fstar(x,1,4.5),xlim=c(0,3), xlab="Morgans", ylab="Density")
title("Crossover distance")
abline(h=1,col="gray")





dens = sapply(seq(0,4,by=0.001), function(x) fstar(x,1,3.2) )

cdf = cumsum(dens)
cdf = cdf/max(cdf)

cdf[1:10]
crossoverCDFvector = cdf


generateRecombinationDistances = function ( n ) {
    
    t = findInterval( runif(n), crossoverCDFvector )
    t = t * 4 / length(crossoverCDFvector)
    100 * t
}

hist(generateRecombinationDistances(20000),n=100)




#### Other - not used ####


getBpFromCm = function(cm) {
    ind = findInterval(cm,geneticMap$cm) 
    if(ind > 10) {
        s = (ind-10):(ind+10)
    } else { 
        s = 1:20
    }
    
    approx(geneticMap$cm[s], geneticMap$bp[s], cm)
}




getCmFromBp = function(bp) {
    ind = findInterval(bp,geneticMap$bp) 
    if(ind > 10) {
        s = (ind-10):(ind+10)
    } else { 
        s = 1:20
    }
    
    approx(geneticMap$bp[s], geneticMap$cm[s], bp)
}



system.time( x <- lapply(seq(1e6,10e6,l=1000), getCmFromBp) )
system.time( x <- approx( geneticMap$bp, geneticMap$cm, seq(1e6,10e6,l=1000000) ) ) 


chromLength = max(a[,4])



range(a$cm)




t=findInterval(runif(100000),cdf)
hist(t,n=100)    

mean(t)



f = function() {
    
    
    
    pos = generateRecombinationDistances(50)
    pos
    pos = cumsum(pos)
    pos
    maxCm = 1200
    pos = pos [ pos <= maxCm ]
    if(length(pos) < 2) pos = c(pos,99999999)
    
    
    recombVector = rep(0, 101)
    recombPositions = seq(0,30,l=length(recombVector))
    
    
    for(i in 1:length(pos)) {
        s = recombPositions > pos[i] 
        #if(runif(1) < 0.5) 
            recombVector [ s ]  = recombVector [ s ]  + 1
    }
    recombVector = recombVector %% 2
    
    recombVector 
    recombVector[length(recombVector)]
    
    
    
    
}




(1-exp(-2*30/100))/2

mean ( replicate(31010,f()) )





