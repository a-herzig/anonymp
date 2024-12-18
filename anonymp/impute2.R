###Notebook of commands to create Figure 4, in ANONYMP paper, detailing a simple HMM
###Anthony Herzig Dec 2024
###anthony.herzig@inserm.fr

###Materials can be found in https://lysine.univ-brest.fr/anonymp/
###R commands

###Read the reference panel with Gaston
library(matrixStats)
library(gaston)
library(gaston.utils)
filename<-paste("/PROJECTS/ANONYMP/chr15_5popSim_4B11_Ref.vcf.gz",sep="")
###45000 individuals, read the column names
namesChromR<-strsplit(system(paste("zcat ",filename," | grep -n -m 1 \"CHROM\"",sep=""),intern=T),":")[[1]][2]
namesChromR<-unlist(strsplit(namesChromR,"\t"));namesChromR[1]<-"CHROM";samples<-namesChromR[-c(1:9)]
###45000 individuals, read the haplotypes
REF_haplo_raw <- .Call("gg_read_vcf_chr_range_haplo", PACKAGE = "gaston.utils", filename, FALSE, -1L, 0L, 0L, samples)
###Read the genotypes, but more importantly the variant data, positions etc.
REF <- read.vcf(file=filename)
chr<-15
#### carte de recombination
map <- read.table(paste("/PUBLIC_DATA/GeneticMaps/Bherer/Refined_EUR_genetic_map_b37/sexavg_chr15_s2.txt",sep=""),header=T,as.is=T)
arrayS<-read.table("/PROJECTS/ANONYMP/chip.txt",header=F,as.is=T)
arrayS[,1]<-15

### interpolation of the recombination distances onto the list of positions in our data
m <- match(REF@snps$pos, map[,1]) 
distcM<-cbind(REF@snps$pos,map[m,3],rep(1,length(REF@snps$pos)))
map2 <- map[is.na(match(map[,1],distcM[,1])),c(1,3)] ; map2 <- cbind(map2,rep(2,nrow(map2)))
c("pos","dist","tag") -> colnames(distcM) -> colnames(map2)
map_tagged <- rbind(distcM,map2)
map_tagged <- map_tagged[order(map_tagged[,1], decreasing = FALSE),]
w2<-which(is.na(map_tagged[,2]))
d3<-approx(as.numeric(map_tagged[,1]),as.numeric(map_tagged[,2]),xout=map_tagged[w2,1],rule=1)
map_tagged[w2,2]<-d3$y
distcM <- map_tagged[which(map_tagged[,3] == 1),c(1,2)]
###store the distances onto the REF object (the reference panel)
REF@snps$dist <- distcM[,2]

m1<-lm(REF@snps$dist~REF@snps$pos)
REF@snps$dist[is.na(REF@snps$dist)]<-m1$fitted.values[is.na(REF@snps$dist)]-max(m1$fitted.values[is.na(REF@snps$dist)])+min(REF@snps$dist,na.rm=T)
plot(REF@snps$dist)

rm(list=c("map2","map_tagged","w2","d3","m"))

###form a gaston object with the haplotypes - gaston thinks I have 90000 individuals, in fact we have 90000 haplotypes
ped <- data.frame(famid = 1, id = paste0(rep(samples,each=2),c("_h1","_h2")), father = 0, mother = 0, sex = NA, pheno = NA, stringsAsFactors = FALSE)
snp <- data.frame(chr = 15, id = ifelse( REF_haplo_raw$id ==".", paste(REF_haplo_raw$chr,REF_haplo_raw$pos,REF_haplo_raw$A1,REF_haplo_raw$A2,sep=":"),REF_haplo_raw$id), dist = REF@snps$dist, pos = REF_haplo_raw$pos , A1 = REF_haplo_raw$A1, A2 = REF_haplo_raw$A2,quality = REF_haplo_raw$quality, filter = factor(REF_haplo_raw$filter), stringsAsFactors = FALSE)
REF_haplo <- new("bed.matrix", bed = REF_haplo_raw$bed, snps = snp, ped = ped, p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )


####Do the same for the Target data: 2000 extra individuals - Cases = Target
filename<-paste("/PROJECTS/ANONYMP/chr15_5popSim_4B11_Cases_FULL.vcf.gz",sep="")
namesChromR<-strsplit(system(paste("zcat ",filename," | grep -n -m 1 \"CHROM\"",sep=""),intern=T),":")[[1]][2]
namesChromR<-unlist(strsplit(namesChromR,"\t"));namesChromR[1]<-"CHROM";samples<-namesChromR[-c(1:9)]
TARGET_haplo_raw <- .Call("gg_read_vcf_chr_range_haplo", PACKAGE = "gaston.utils", filename, FALSE, -1L, 0L, 0L, samples)


pedT <- data.frame(famid = 1, id = paste0(rep(samples,each=2),c("_h1","_h2")), father = 0, mother = 0, sex = NA, pheno = NA, stringsAsFactors = FALSE)
snpT <- snp[which(REF_haplo@snps$pos%in%arrayS[,2]),]

TARGET_haplo <- new("bed.matrix", bed = TARGET_haplo_raw$bed, snps = snp, ped = pedT, p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE, standardize_mu_sigma = FALSE )


refHap<-as.matrix(REF_haplo)
target<-as.matrix(TARGET_haplo)
### note I will use rownames and colnames of these objects later

w1<-which(REF_haplo@snps$pos%in%arrayS[,2])


###form my haplotype matrices on the positions in chip.txt
H<-t(refHap[,w1])
G<-t(target[,w1])


###how to split the chromosme into 16 'chunks'
chunks<-read.table("/PROJECTS/ANONYMP/coordinates.I51.2_15.txt",header=F,as.is=T)

### select the best K haplotypesn from matrix H for imputing haplotye i in matrix G
Hamm<-function(i,H,G,K){
dh1<-colSums(H==G[,i])
return(order(dh1,decreasing=T)[1:K])
}

###same-diff matrix for haplotype i in matrix G and their closest K haplotypes in matrix H
accord<-function(i,H,G,K){
dh1<-Hamm(i,H,G,K)
mh1<-1*(H[,dh1]==G[,i])
colnames(mh1)<-dh1
return(mh1)
}

###Set K here
K<-500
rho<- 1-exp((-40000)*diff(REF@snps$dist[w1]/100)/K)
tmat<-log(cbind((1-rho) + rho/K,rho/K))

emit<-c(log(0.0001),log(0.9999))


###chunkwhich to select which snps to use for a chunk and stock the haplotypes needed and the distances needed
chunkwhich<-function(i,chunks){

s_e<-unlist(strsplit(chunks[i,4],":"))[2]
s_e<-as.numeric(unlist(strsplit(s_e,"-")))

s_e2<-unlist(strsplit(chunks[i,3],":"))[2]
s_e2<-as.numeric(unlist(strsplit(s_e2,"-")))


w<-which(snpT$pos>= s_e[1] & snpT$pos<= s_e[2])
w<-c((min(w)-250):(min(w)-1),w,(max(w)+1):(max(w)+250))
w<-w[which(w>=1)];w<-w[which(w<=length(snpT$pos))]
w2<-which(REF_haplo@snps$pos>= s_e2[1] & REF_haplo@snps$pos<= s_e2[2])

l1<-list(snpT$id[w],REF_haplo@snps$id[w2],H[w,],G[w,],length(w),tmat[w[1:(length(w)-1)],])
names(l1)<-c("snpids","refsnpids","Hc","Gc","nSNP","tmatc")
return(l1)
}


forward_backward<-function(i,lcj,K){

mh1<-accord(i,lcj$Hc,lcj$Gc,K)
N<-lcj$nSNP
alphaM <- matrix(NA, ncol = N, nrow = K)
alphaM[,1] <- (-log(K)) + emit[mh1[1, ]+1]
alphaM[,1]<-alphaM[,1]-logSumExp(alphaM[,1])
betaM <- matrix(NA, ncol = N, nrow = K)
betaM[,N]<-0
betaM[,N]<-betaM[,N]-logSumExp(betaM[,N])

for(j in 2:N){
###log sum of the previous alpha values execpt for one value
abar<-sapply(1:K,function(i) logSumExp(alphaM[-i,j-1]))
###x1 is the probabilities of staying in the same state 
x1<-lcj$tmatc[j-1,1]+alphaM[,j-1]
###x2 is the proabilities of moving to any other state
x2<-lcj$tmatc[j-1,2]+abar
alphaM[,j]<-sapply(1:K,function(i) logSumExp(c(x1[i],x2[i])))+emit[mh1[j,]+1]
alphaM[,j]<-alphaM[,j]-logSumExp(alphaM[,j])

b1 <- emit[mh1[N-j+1,]+1] + betaM[,N-j+2]
bbar<-sapply(1:K,function(i) logSumExp(b1[-i]))
x1<-lcj$tmatc[j-1,1]+b1
x2<-lcj$tmatc[j-1,2]+bbar
betaM[,lcj$nSNP-j+1]<-sapply(1:K,function(i) logSumExp(c(x1[i],x2[i])))
betaM[,lcj$nSNP-j+1]<-betaM[,lcj$nSNP-j+1]-logSumExp(betaM[,lcj$nSNP-j+1])
}



logppm <- alphaM + betaM
logppm <- sweep(logppm, 2, colLogSumExps(logppm), "-")
rownames(logppm)<-REF_haplo@ped$id[as.numeric(colnames(mh1))]
colnames(logppm)<-lcj$snpids

allsnps<-as.character(sort(as.numeric(unique(c(lcj$snpids,lcj$refsnpids)))))
wA<-which(allsnps%in%lcj$refsnpids)
wA<-c(min(wA)-1,wA,max(wA)+1)
wA<-wA[which(wA>=1)];wA<-wA[which(wA<=length(allsnps))]

logppmF<-matrix(NA,nrow=nrow(logppm),ncol=length(wA))
rownames(logppmF)<-rownames(logppm)
colnames(logppmF)<-allsnps[wA]
logppmF[,colnames(logppm)[which(colnames(logppm)%in%colnames(logppmF))]]<-logppm[,which(colnames(logppm)%in%colnames(logppmF))]

return(as.matrix(logppmF))

}



###I5 is impute5 results, M4 is minimac4 results
I5<-read.table("/PROJECTS/ANONYMP/SAP.I5_anonympinfo.temp",header=F,as.is=T)
I5<-cbind(I5,read.table("/PROJECTS/ANONYMP/SAP.I5_anonymp.temp",header=F,as.is=T))
rownames(I5)<-REF_haplo@snps$id;I5<-I5[,-c(1:9)]
colnames(I5)<-paste0(rep(samples,each=2),c("_h1","_h2"))

M4<-read.table("/PROJECTS/ANONYMP/SAP.M4_anonympinfo.temp",header=F,as.is=T)
M4<-cbind(M4,read.table("/PROJECTS/ANONYMP/SAP.M4_anonymp.temp",header=F,as.is=T))
rownames(M4)<-REF_haplo@snps$id;M4<-M4[,-c(1:9)]
colnames(M4)<-paste0(rep(samples,each=2),c("_h1","_h2"))


DT<-c()
TT<-c()
IT<-c()
MT<-c()

linear_interpolation<-function(x){
return(approx(x=REF_haplo@snps$dist[which(colnames(refHap)%in%colnames(ppm))],y=as.numeric(exp(x)),xout=REF_haplo@snps$dist[which(colnames(refHap)%in%colnames(ppm)[wk])],rule=2)$y)
}


###in practice this will be massivly parallelised
for(chunk in 1:16){
lcj<-chunkwhich(chunk,chunks)
print(chunk)
dc<-c()
tc<-c()
ic<-c()
mc<-c()

for(z in 1:4000){
print(z)
ppm<-forward_backward(z,lcj,K)
wk<-which(is.na(ppm[1,]))

m1<-apply(ppm,1,linear_interpolation)

dc<-cbind(dc,colSums(t(m1)*refHap[rownames(ppm),colnames(ppm)[wk]]))
tc<-cbind(tc,target[z,colnames(ppm)[wk]])
ic<-cbind(ic,I5[colnames(ppm)[wk],colnames(G)[z]])
mc<-cbind(mc,M4[colnames(ppm)[wk],colnames(G)[z]])
}

DT<-rbind(DT,dc)
TT<-rbind(TT,tc)
IT<-rbind(IT,ic)
MT<-rbind(MT,mc)
}





save(TT,file="/PROJECTS/ANONYMP/TT.Rdata")
save(DT,file=paste0("/PROJECTS/ANONYMP/DT_k",K,".Rdata"))
save(IT,file="/PROJECTS/ANONYMP/IT.Rdata")
save(MT,file="/PROJECTS/ANONYMP/MT.Rdata")




#####

load(file="/PROJECTS/ANONYMP/TT.Rdata")
load(file="/PROJECTS/ANONYMP/DT_k100.Rdata")
load(file="/PROJECTS/ANONYMP/DT_k500.Rdata")
load(file="/PROJECTS/ANONYMP/IT.Rdata")
load(file="/PROJECTS/ANONYMP/MT.Rdata")



###Make Figure 4

RT<-sapply(1:ncol(TT),function(x) cor(DT[,x],TT[,x])**2)
RT<-cbind(RT,sapply(1:ncol(TT),function(x) cor(DT500[,x],TT[,x])**2))
RT<-cbind(RT,sapply(1:ncol(TT),function(x) cor(IT[,x],TT[,x])**2))
RT<-cbind(RT,sapply(1:ncol(TT),function(x) cor(MT[,x],TT[,x])**2))


par(mfrow=c(1,3))
plot(RT[,1:2],xlim=c(0.89,1),ylim=c(0.94,1),xlab="R2 ANONYMP K=100",ylab="ANONYMP K=500",main="(a) R2 per haplotype");abline(0,1,col="red")
plot(RT[,c(2,3)],xlim=c(0.89,1),ylim=c(0.94,1),xlab="R2 ANONYMP K=500",ylab="R2 IMPUTE5",main="(b) R2 per haplotype");abline(0,1,col="red")
plot(RT[,c(2,4)],xlim=c(0.89,1),ylim=c(0.94,1),xlab="R2 ANONYMP K=500",ylab="R2 MINIMAC4",main="(c) R2 per haplotype");abline(0,1,col="red")



RT2<-cbind(sapply(1:nrow(TT),function(x) cor(DT[x,],TT[x,])**2),sapply(1:nrow(TT),function(x) cor(DT500[x,],TT[x,])**2),sapply(1:nrow(TT),function(x) cor(IT[x,],TT[x,])**2),sapply(1:nrow(TT),function(x) cor(MT[x,],TT[x,])**2))

maf1<-rowMeans(TT)

l1<-c(0,0.005,0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.99,1)
RT2b<-c()
lan1<-c()
for(j in 2:length(l1)){
lan1<-c(lan1,paste0("(",l1[j-1],"-",l1[j],"]"))
w<-which(maf1<=l1[j]&maf1>l1[j-1])
RT2b<-rbind(RT2b,colMeans(RT2[w,],na.rm=T))
}


par(mfrow=c(1,1))
plot(RT2b[,1],ylim=c(0.78,1),typ="b",xaxt="n",xlab="True Minor Allele Frequency",ylab="aggregate R2",pch=15,lty=2,cex=1.5,main="(d) aggregate R2 across MAF bins")
points(RT2b[,2],typ="b",pch=17,lty=3,cex=1.5)
points(RT2b[,3],typ="b",pch=18,lty=4,cex=1.75)
points(RT2b[,4],typ="b",pch=19,lty=5,cex=1.5)
legend(x=3.5,y=0.84,c("ANONYMP K=100","ANONYMP K=500","IMPUTE5","MINIMAC4"),pch=c(15,17,18,19),lty=2:5,pt.cex=c(1.5,1.5,1.75,1.5))
axis(1,at=1:10,labels=lan1)


