library(hapsim)

data(ACEdata)

# create the haplotype object 
str(ACEdata)


t = ACEdata
image(as.matrix(cbind(ACEdata)))


npop = 100
t = data.frame(s1 = rep(0,npop) )
f1 = runif(npop,0,0.01)
for(i in 1:npop)
t[,paste(i)] = rbinom(npop,1,f1[i])

t = t [ , apply(t,2,max) > 0 ]

x <- haplodata(t) 

# simulates a first sample of 100 haplotypes using all markers
y1 <- haplosim(9000, x) 
# compares allele frequencies in real and simulated samples
plot(x$freqs, y1$freqs, title=paste("MSE:",y1$mse.freqs)); abline(a=0, b=1)

# compares LD coefficients in real and simulated samples
ldplot(mergemats(x$cor, y1$cor), ld.type='r') 

# simulates a second sample of 1000 haplotypes using the first 20 markers only
y2 <- haplosim(1000, which.snp=seq(30), x) 


