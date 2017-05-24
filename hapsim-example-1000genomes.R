library(hapsim)

data(ACEdata)

# create the haplotype object 
str(ACEdata)

range(ACEdata)
t = ACEdata
image(as.matrix(cbind(ACEdata)))


npop = 100
nmark = 1500

t = data.frame(s1 = rep(0,npop) )

f1 = runif(nmark,0,0.41)

for(i in 1:nmark)
t[,paste(i)] = rbinom(npop,1,f1[i])



t[,10:20] = t[,30:40]
t = t [ , apply(t,2,max) > 0 ]


t = t(cbind( gt1[1:400,], gt2[1:400,]))

x <- haplodata(t) 


# simulates a first sample of 100 haplotypes using all markers
y1 <- haplosim(1000, x) 


q = y1$data

haplo1 = apply(q,1,function(x) paste(x,collapse=""))
haplo2 = apply(t,1,function(x) paste(x,collapse=""))


haplo1 %in% haplo2



str(y1)

# compares allele frequencies in real and simulated samples
plot(x$freqs, y1$freqs, title=paste("MSE:",y1$mse.freqs)); abline(a=0, b=1)

# compares LD coefficients in real and simulated samples
ldplot(mergemats(x$cor, y1$cor), ld.type='r') 




x <- haplodata(t) 

t2 = as.matrix(t)

i = sample(1:22,311,replace=T)
t3 = t2[i,]
x2 = haplodata(t3)
str(x2)


ldplot(mergemats(x$cor, x2$cor), ld.type='r') 














# simulates a second sample of 1000 haplotypes using the first 20 markers only
y2 <- haplosim(1000, which.snp=seq(30), x) 
