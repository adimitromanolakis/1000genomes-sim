

#  sh R-1000genomes-extract.sh 1000000-1400000 1 area.vcf CEU.txt

rm(list=ls())

vcf = read.table("area.vcf",sep="\t",comment="#",as=T,h=T)
str(vcf)
vcf = vcf[,-ncol(vcf)]

vcf = vcf[1:100,]
str(vcf)


BP = vcf[,2]

plot(diff(BP),t="l")

gt = vcf[,-(1:3)]
gt[1:5,1:5]
#table(unlist(gt[,]))


gt = as.matrix(gt)
gt = gt[,-ncol(gt)]
gt[is.na(gt)] = "0|0"



# Remove multiallelic sites

goodGT = c("0|0","0|1","1|0","1|1")
badSites = apply(gt,1,function(x) sum(!(x %in% goodGT)))

gt = gt[badSites == 0,]
      
table(gt)


table( gt[3,] =="0|0" ,useNA="a")





maf = apply(gt,1, function(x) mean(x != "0|0"))

gt = gt[ maf>0.01 & maf < 0.99,]
#gt = gt[ maf>0.1 & maf < 0.9,]



maf = apply(gt,1, function(x) mean(x != "0|0"))
range(maf)

dim(gt)

str(gt)



library(stringr)


haplo = apply(gt,2, function(x) paste(x,sep="",collapse=" "))
data.frame( table(haplo) )

gt1 = t(apply(gt,1,function(x) str_sub(x,1,1)))
gt2 = t(apply(gt,1,function(x) str_sub(x,3,3)))

gt1 = t(apply(gt1,1,as.numeric))
gt2 = t(apply(gt2,1,as.numeric))


haplo1 = apply(gt1,2, function(x) paste(x,sep="",collapse=" "))
haplo2 = apply(gt2,2, function(x) paste(x,sep="",collapse=" "))

haplos = rbind(haplo1, haplo2)


dim(gt)
dim(gt1)
dim(gt2)

str(gt1)

library(gplots)

m = apply(gt1,1,as.numeric)
m = m + apply(gt2,1,as.numeric)
#m = m[,1:600]



dim(m)
h1 = function(x) hclust(x,method="average")
d1 = function(x) dist(x, method="manh")

heatmap.2(m,   Rowv=T,Colv=T,trace="n",hclustfun = h1,distfun=d1)









colorPal = function(pal = "RdYlBu", n=300) {
    
    library(RColorBrewer)
    
    #  c1 = colorRampPalette(brewer.pal(300,"Blues"))
    #  c1 = colorRampPalette(brewer.pal(300,"YlOrRd"))
    c1 = colorRampPalette(brewer.pal(n,pal))
    
    c1(n)
    
}


cor1 = cor(m)
c1 = colorRampPalette(c("#ffffff",brewer.pal(99,"YlOrRd")))

heatmap.2(cor1[1:500,1:500]^2,trace="n",col=c1(200), Rowv=F,Colv=F)

dim(gt)








orderRow = (h1(d1(m))$order)
orderCol = rev(h1(d1(t(m) ))$order)
str(orderCol)
dim(m)




heatmap.2(m[orderRow,orderCol],    Rowv=T,Colv=F,trace="n",hclustfun = h1,distfun=d1)

m = apply(gt2,1,as.numeric)
#m = m + apply(gt2,1,as.numeric)
m = m[,1:600]
heatmap.2(m[orderRow,orderCol],    Rowv=T,Colv=F,trace="n",hclustfun = h1,distfun=d1)






haplo = apply(gt1,2, function(x) paste(x[1:66],sep="",collapse=" "))
data.frame( table(haplo) )

haplo = apply(gt2,2, function(x) paste(x[1:66],sep="",collapse=" "))
data.frame( table(haplo) )



