library(hapsim)
library(compiler)
enableJIT(3)


maxNumberOfVariants = 400
min_maf = 0.25

randomdata = 0



# node ../fq-thin-vcf-file.js --out vt.vcf --in chr2.vcf --thin 0.7

# bcftools view -s NA19904,NA19913,NA20317,NA20318,NA20334,NA20355,NA20359,NA20362 vt.vcf|grep -v '^##' > vt2.vcf

#4 90645249 90759446 strand -1 id SNCA name SNCA bioType protein_coding
'
CHROM=4
VCF=~/1000genomes/ALL.chr"$CHROM".*.vcf.gz
REGION=4:90645249-90759446
REGION=4:60645249-61569446

bcftools view -S CEU.txt --force-samples -r $REGION  $VCF > /tmp/1.vcf
grep -v "^##" /tmp/1.vcf > haplosims/1.vcf



'

if(class(vcf) != "data.frame") {
    
#### Read VCF ####

library(stringr)


cat("[#.......] Reading VCF file..\n")


vcf=1


vcf = read.table("haplosims/1.vcf",sep="\t",comment="",as=T,h=T)

vcf = vcf[str_length(vcf[,5]) == 1 ,]


# Reduce number of variants

vcf = vcf[seq(1,nrow(vcf),by=4) ,]


cat("[##......] Chromosome:  ", unique(vcf[,1]), 
    " Mbp: " , min(vcf[,2])/1e6, 
    " Region Size: ", 1e-3 * ( max(vcf[,2])-min(vcf[,2]) ) ,"kb ",
    "Num of variants:", dim(vcf)[1] ,
    "\n");


#vcf = vcf[10000:10500,]


library(stringr)

gt = vcf[,-(1:9) ]

#gt[1:5,1:5]


gt1 = apply( gt, 2,  function(x) as.numeric( str_sub(x,1,1) )  )
gt2 = apply( gt, 2,  function(x) as.numeric( str_sub(x,3,3) )  )


#dim(gt)
#dim(gt1)
#dim(gt2)



#### --- Filter by MAF ---- ####


cat("[###.....] Filtering and thinning variants\n");

maf = apply(gt1,1,function(x) mean(x,na.rm=T))
maf[maf>0.5] = 1 - maf[maf>0.5]


maf2 = apply(gt2,1,function(x) mean(x,na.rm=T))
maf2[maf2>0.5] = 1 - maf2[maf2>0.5]

ok = apply(gt1,1,function(x) max(x,na.rm=T))
s = which(maf > min_maf & maf2 > min_maf & ok < 2)


if(length(s)>maxNumberOfVariants) s = sort( sample(s,maxNumberOfVariants) )


gt1 = gt1[s,]
gt2 = gt2[s,]
gt = gt[s,]
vcf = vcf[s,]

# 


cat("## Chromosome:  ", unique(vcf[,1]), 
    " Mbp: " , min(vcf[,2])/1e6, 
    " Region Size: ", 1e-3 * ( max(vcf[,2])-min(vcf[,2]) ) ,"kb ",
    "Num of variants:", dim(vcf)[1] ,
    "\n");

length(s)


dim(gt1)
dim(vcf)


#s = sort( sample(s) )



##

cat("[####....] Reading genetic map\n");

readGeneticMap = function(chromosome) {
    
        filelocation = sprintf( "geneticmap/genetic_map_GRCh37_chr%s.txt", as.character( chromosome) ) 
        a = read.table(filelocation, as=TRUE, h=TRUE)
        colnames(a) = c("chr","bp","rate.cm","cm")
        a = a[ order(a$cm),]
        a
        
}

geneticMap = readGeneticMap(4)


##

bp = vcf[,2]
cm = approx( geneticMap$bp, geneticMap$cm, bp )$y

plot(bp/1e6,cm, t="l",cex=0.2)
segments( bp/1e6,cm,  bp/1e6,cm -0.2, lwd=0.3,col="blue")



}

#### Mating and recombination #####




generateRecombinationDistances = function ( n ) {
    
    ( rexp(n, 0.01))
}




# Genetates a vector of origin of segments (0 - haplotype1 ,  1 - haplotype2 )

generateSingleRecombinationVector = function(cm) {
    
    n = length(cm)
    maxCM = max(cm)
    
    pos = generateRecombinationDistances(14)
    
    while(sum(pos) < maxCM) pos = c(pos, generateRecombinationDistances(10) )
    if(sum(pos) < maxCM) { stop("Not enough recombination events"); }
    
    pos = cumsum(pos)
    
    recombVector = rep(TRUE, n)
    
    p = findInterval(pos,cm)
    
    
    for(i in 1:length(pos)) {
        
        if(p[i] < n) recombVector[ p[i]:n ] =  !recombVector[ p[i]:n ]  
        
        #s = cm > pos[i] 
        #if(runif(1) < 0.5) 
        #recombVector [ s ]  = recombVector [ s ]  + 1
        
    }
    
    #  recombVector = recombVector %% 2
    if(runif(1,0,1) < 0.5) recombVector = ! recombVector
    
    recombVector+0
    
}

#generateSingleRecombinationVector(cm*400)





#hist(generateRecombinationDistances(130000),n=100)

if(0) {
recombVector = generateSingleRecombinationVector(cm)
paste(recombVector[seq(1,length(recombVector),by=50)],collapse="")




n = length(cm)
haplo3 = rep(0, n)




cm = seq(0,331,l = length(1:1000) )

ibd2 = function() {
    
    origin3father = generateSingleRecombinationVector(cm)
    
    origin3mother =  generateSingleRecombinationVector(cm)
    
    origin4father = generateSingleRecombinationVector(cm)
    
    origin4mother = generateSingleRecombinationVector(cm)
    
    mean(origin3father == origin4father | origin3mother == origin4mother) 
    
}



system.time ( t <- replicate(300,ibd2()) )
mean(t)
}

generateSingleRecombinationVector = cmpfun(generateSingleRecombinationVector)




# PLOT: Recombination points for 200 individuals across region
if(0) {
        cm = seq(0,111,l = length(1:1000) )
        
        n1=length(cm)
        n2 = 200
        m = matrix(NA,nrow=n1,ncol=n2)
        m=apply(m,2,function(x) generateSingleRecombinationVector(cm))
        
        #m[,seq(1,n2,by=2)  ] = 0
        
        image(m)
}





# PLOT: Mean number of crossover points in a single region per generation
#hist(replicate(11900, sum(diff(generateSingleRecombinationVector(seq(0,300,l=120)) ) != 0  ) ), breaks=0:12)

# SYSTEM: Timing
#system.time(replicate(1000,   generateRecombinationDistances(10)       ))
#system.time(replicate(1000,   generateSingleRecombinationVector(cm)    ))


if(0) {
# Compute mean number of recombinations per generation

cm = seq(0,25,l=200)
v = generateSingleRecombinationVector(cm)
paste(v,collapse="")

mean( sapply(1:5000, function(i) { v = generateSingleRecombinationVector(cm); abs( max(v)-min(v) ); } ) )
}



#### Haplosim 2 function ####


haplosim2 = function (n, hap, which.snp = NULL, seed = NULL, force.polym = TRUE, 
                      summary = TRUE) 
{
    if (length(seed) > 0) 
        set.seed(seed)
    nsubset <- length(which.snp)
    if (nsubset > 0) {
        which.snp <- sort(unique(which.snp))
        if ((which.snp[1] < 1) || (which.snp[nsubset] > length(hap$freqs))) 
            stop("which.snp does not contain valid indeces")
        indexset <- which.snp
    } else { indexset <- seq(length(hap$freqs)) }
    
    nloci <- length(indexset)
    quants <- qnorm(hap$freqs[indexset])
    y <- matrix(0, nrow = n, ncol = nloci)
    colnames(y) <- names(hap$freqs[indexset])
    A <- mvrnorm(n, mu = rep(0, nloci), Sigma = hap$cov[indexset, 
                                                        indexset], tol = 1e-06, empirical = FALSE)
    if (n == 1) 
        A <- array(A, dim = c(1, nloci))
    for (i in 1:n) {
        for (j in 1:nloci) {
            if (A[i, j] <= quants[j]) 
                y[i, j] <- 0
            else y[i, j] <- 1
        }
    }
    y.freqs <- allelefreqs(y)
    
    y.freqs <- y.freqs$freqs
    if (summary) {
        y.cor <- cor(y)
        y.div <- divlocus(y)
        mse.freqs <- mse(hap$freqs[indexset], y.freqs)
        hap.cor <- hap$cor[indexset, indexset]
        tmat1 <- hap.cor[upper.tri(hap.cor)]
        tmat2 <- y.cor[upper.tri(y.cor)]
        mse.cor <- mse(tmat1, tmat2)
        return(list(data = y, freqs = y.freqs, cor = y.cor, div = y.div, 
                    mse.freqs = mse.freqs, mse.cor = mse.cor))
    }
    else return(list(data = y, freqs = y.freqs))
}




##### Fake Geno


######## SIM data object #####

cat("[#####...] Creating SIM object\n");

SIM = new.env()


# Original haplotypes from 1000 genomes / other data

SIM$population_gt1 = gt1
SIM$population_gt2 = gt2



if(randomdata) {
    
    str(SIM$population_gt1)
    # indoviduals in columns
    
    x = SIM$population_gt1*0
    for(i in 1:nrow(x)) x[i,] = rbinom( ncol(x), 1 , 0.25)
    SIM$population_gt1 = x
    SIM$population_gt2 = x
    
    SIM$haplodata = haplodata(  t( cbind(SIM$population_gt1, SIM$population_gt2) ) ) 
    
#    ldplot(SIM$haplodata$cor,ld.type="r")
    
    
    
}




# Generate haplodata object 

dim( t( cbind(SIM$population_gt1, SIM$population_gt2) )  )
SIM$haplodata = haplodata(  t( cbind(SIM$population_gt1, SIM$population_gt2) ) ) 

#SIM$N_people = nrow(SIM$population_gt1)
#SIM$N_variants = ncol(SIM$population_gt1)

# Variant information

SIM$varinfo = vcf[,1:8]
SIM$bp = vcf[,2]

SIM$cm = approx( geneticMap$bp, geneticMap$cm, SIM$bp )$y


SIM$N_markers = nrow(gt1)

SIM$pool = NA
SIM$npool = 0



SIM$total_individuals = 400
SIM$individuals_generated = 0


SIM$gt1 = matrix(nrow = SIM$total_individuals, ncol = SIM$N_markers, -1)
SIM$gt2 = matrix(nrow = SIM$total_individuals, ncol = SIM$N_markers, -1)






SIM$generateNewHaplotypes = function(n) {
    
    if(SIM$npool < 2) {
        cat("Generate new individual pool n=100\n");
        
        SIM$pool <<- haplosim2(200, SIM$haplodata, summary = F)$data
        SIM$npool <<- 200
        
    }
   
   GT = list(gt1 =  SIM$pool[SIM$npool,], gt2 = SIM$pool[SIM$npool-1,] )
   SIM$npool <<- SIM$npool - 2
   
   return(GT)
   
}


addUnrelatedIndividual = function() {

    newGenotypes = SIM$generateNewHaplotypes()
    
    if(SIM$individuals_generated >= SIM$total_individuals) {
        
     stop("No more space for saving new individual genotypes")
        
    }
    
    
    SIM$individuals_generated = SIM$individuals_generated + 1
    j = SIM$individuals_generated
    
    cat("Adding individual ",j, " from pool\n");
    
    
    SIM$gt1[j,] = newGenotypes$gt1
    SIM$gt2[j,] = newGenotypes$gt2
    
    return(j)
}


addIndividualFromGenotypes = function(gt1,gt2) {
    
    if(SIM$individuals_generated >= SIM$total_individuals) {
        
        stop("No more space for generating new individual genotypes")
        
    }

    SIM$individuals_generated = SIM$individuals_generated + 1
    j = SIM$individuals_generated
    
    cat("Adding individual ",j, " from specified genotypes\n");
    
    SIM$gt1[j,] = gt1
    SIM$gt2[j,] = gt2
    
    return(j)
}

newNuclearFamily = function(fid) {
    
    fam = data.frame(fid = fid  , 
                     id = c(1,2,3) , 
                     father = c(0,0,1) , 
                     mother = c(0,0,2) , 
                     sex = c(1,2,1)
                )
    
    
    j1 = addUnrelatedIndividual()
    j2 = addUnrelatedIndividual()
    j3 = mate(j1,j2)
    
    
    fam$gtindex = c(j1,j2,j3)
    
    fam
}


newFamilyWithOffspring = function(fid, n = 2) {
    
    fam = data.frame(fid = fid  , 
                     id = c(1:2) , 
                     father = c(0,0), 
                     mother = c(0,0), 
                     sex = c(1,2)
    )
    
    
    j1 = addUnrelatedIndividual()
    j2 = addUnrelatedIndividual()
    
    fam$gtindex = c(j1,j2)
    
    
    
    
    
    for(i in 1:n) {
     j3 = mate(j1,j2)
     newFamLine = c(fid, i+10, 1,2, 1 , j3)
     fam = rbind(fam, newFamLine)
    }
    
    return (fam)
}


#SIM$cm = seq(0, 3200, l=dim(SIM$gt1)[2] )
SIM$npool = 0

mate = function(i,j) {
    
    # For now, we hope i,j are of different sex
    
    recomb1 = generateSingleRecombinationVector(SIM$cm)
    recomb2 = generateSingleRecombinationVector(SIM$cm)
    
    
    # FATHER1 = SIM$gt1[i,]
    # FATHER2 = SIM$gt2[i,]
    # MOTHER1 = SIM$gt1[i,]
    # MOTHER2 = SIM$gt2[i,]
    
    
    
    gt1 =  recomb1 * SIM$gt1[i,]  +  (1-recomb1) * SIM$gt2[i,]
    gt2 =  recomb2 * SIM$gt1[j,]  +  (1-recomb2) * SIM$gt2[j,]
    
    #gt1 = SIM$gt1[i,]
    #gt1[recomb1==1] = SIM$gt2[i,][recomb1==1]
    
    #gt2 = SIM$gt1[j,]
    #gt2[recomb1==1] = SIM$gt2[j,][recomb1==1]
 
    #cat("Mate: ",i," " ,j,"\n");
    
    #p1 = paste(SIM$gt1[i,],collapse="")
    #p2 = paste(SIM$gt2[i,],collapse="")
    #p3 = paste(gt1,collapse="")
    
    #cat("R1:", paste(recomb1,collapse=""),"\n",sep="")
    #cat("R2:", paste(recomb2,collapse=""),"\n",sep="")
    
    #cat(p1,"\n",p2,"\n",p3,"\n",sep="")

    if(runif(1) < 0.5) { t = gt1; gt1 = gt2; gt2 = t; }

    
    index = addIndividualFromGenotypes(gt1, gt2)
    
    last_gt <<- gt1 + gt2

    return (index)    

}




pIBS0 = function(i,j) {
    
    x1 = SIM$gt1[i,]
    x2 = SIM$gt2[i,]
    y1 = SIM$gt1[j,]
    y2 = SIM$gt2[j,]
    
    
    v = 2- abs( (x1+x2)- (y1+y2) )
    plot(v,t="l")
    
    table(v)
    sum(v==0)
    
    
}




test1_timingOfMatingFunction = function() {
    
    SIM$individuals_generated <- 0
    addUnrelatedIndividual()
    addUnrelatedIndividual()
    #addUnrelatedIndividual()
    #addUnrelatedIndividual()
    
    
    
    #mate(1,2)
    #g1 = last_gt
    
    for(i in 1:1) {
        print( system.time(x <- sapply(1:50, function(x) mate(1,2))) )
    }
    
}; test1_timingOfMatingFunction()

 
test2_pIBS0_problem = function() {
        SIM$individuals_generated <- 0
    
        ibs0 = c()
        for(i in 1:20) {
        x = newFamilyWithOffspring(2,5)

        
        ind1 = x$gtindex[3]
        ind2 = x$gtindex[5]
        
        v = pIBS0(ind1,ind2)
        ibs0 = c(ibs0, v)
        
        }
        
        cat("IBS0 of full-siblings: " , ibs0, "\n");

}; test2_pIBS0_problem()






#pIBS0(3,4)


# Generate pedigree of 10 individuals

time100families = function() {
    
    SIM$individuals_generated <<- 0
    
    system.time( fam <<- do.call(rbind, lapply(1:25, function(x) newFamilyWithOffspring(x,4) ) ) )
    fam

}

time100families();




str(fam)
dim(SIM$gt1)


s = SIM$haplodata$freqs > 0.2 &  SIM$haplodata$freqs < 0.9
X = (SIM$gt2+SIM$gt1)[1:50,s]






dim(SIM$gt1)
length(SIM$haplodata$freqs)


dim(SIM$gt1)
str(SIM$gt1)

image(abs(cor(t((SIM$gt2+SIM$gt1)[1:50,])))^1.9, col=rev(heat.colors(200)),breaks=seq(0.2,1,l=201))

image(abs(cor( (t(SIM$population_gt1 + SIM$population_gt2)[,])))^1.3, col=rev(heat.colors(200)),breaks=seq(0.2,1,l=201))
image(abs(cor( ((SIM$gt2+SIM$gt1)[1:90,])))^1.3, col=rev(heat.colors(200)),breaks=seq(0.2,1,l=201))



q=(SIM$gt1+SIM$gt2)
table( 2-abs( q[1,]-q[3,]) )


cat("Total individuals generated:", SIM$individuals_generated,"\n");
cat("Genotypes sizes:", dim(SIM$gt1),"\n");





dim(SIM$gt1)
dim(SIM$pool)
dim(SIM$population_gt1)




#### Write ped and map files  ####

N = SIM$individuals_generated 
fam[1,]

simulated_allele1 = SIM$gt1[fam$gtindex ,]
simulated_allele2 = SIM$gt2[fam$gtindex ,]



getGenotypeFromHaplo = function(i) {
    paste( simulated_allele1[i,]+1, simulated_allele2[i,]+1, sep=" ")
}

genotypes = sapply(1:N, getGenotypeFromHaplo)


ids = sprintf("IND%s", seq(1,N))
fam_header = data.frame( fam[,1:4], fam[,5],"1\t",  t(genotypes)) 
write.table(fam_header,file="out.ped",sep=" ",row=F,col=F,quote=F)

dim(fam)
dim(vcf)

map = data.frame(vcf[,1],vcf[,3],0,vcf[,2])
write.table(map,file="out.map",sep=" ",row=F,col=F,quote=F)





system("~/tools/prest-plus/prest --geno out.ped --map out.map --wped --ibs1 0.5")

finished();





#######  Generate random haplotypes #######


if(0) {
r=abs(cor(t(q),use="pair"))
r = r[ nrow(r):1 ,]
heatmap.2(r,trace="n",col=rev(heat.colors(100)),Rowv=F,Colv=F,margins=c(4.2,2),density.info="none",key=F,keysize=0.43)


haplosim = apply(q,1 , function(x) paste(x,collapse=""))
data.frame( table(haplosim) )
}



#####







##### Haplosim will generate simulated individuals in the region ######


if(0) { 
library(hapsim)


dim(gt1)
dim(t(gt1))
z = t(gt1)[,200:900]
z = z[ , apply(z,2,sum) > 0 ]
z = z[ , apply(z,2,mean) < 1  ]



#s = which(apply(z,2,mean) > 0.5)
s
#for(i in s) z[ , i ] = 1- z[,i]

plot( apply(z,2,mean)  ,ylim=c(0,1),cex=0.4)
z = z[,apply(z,2,mean) > 0.1 & apply(z,2,mean) < 0.9] 
      plot( apply(z,2,mean)  ,ylim=c(0,1),cex=0.4)

dim(z)



marker_haplodata = haplodata(z)
str(marker_haplodata)


plot(marker_haplodata$freqs)


s = 1:min(500, ncol(z) )
image(abs(marker_haplodata$cor[s,s]), col=rev(heat.colors(100)),breaks = seq(0,1,l=101))


system.time ( q <- haplosim2(1111,marker_haplodata, summary=T,force.polym = F) )

image(abs(q$cor)[s,s], col=rev(heat.colors(100)),breaks = seq(0,1,l=101))

plot(q$freqs,marker_haplodata$freqs,t="p",pch=20,cex=0.5)



save(gt1,file="/tmp/1.rdata")






image(abs(cor(y)),col=rev(heat.colors(200)) ,breaks=seq(0,1,l=201))
image(abs(marker_haplodata$cor), col=rev(heat.colors(100)),breaks = seq(0,1,l=101))


str( marker_haplodata )
}


