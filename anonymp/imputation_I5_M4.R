###Notebook of commands to perfrom imputation with IMPUTE5 and MINIMAC4 for ANONYMP simulated data
###Anthony Herzig Dec 2024
###anthony.herzig@inserm.fr

###Materials can be found in https://lysine.univ-brest.fr/anonymp/
###R commands

chr<-15

system(paste("bcftools convert -O b -o /PROJECTS/ANONYMP/chr15_5popSim_4B11_Ref.bcf /PROJECTS/ANONYMP/chr15_5popSim_4B11_Ref.vcf.gz",sep=""))
system(paste("tabix -f /PROJECTS/ANONYMP/chr15_5popSim_4B11_Ref.bcf",sep=""))

system(paste("bcftools +fill-tags /PROJECTS/ANONYMP/chr15_5popSim_4B11_Ref.bcf -Ob -o /PROJECTS/ANONYMP/temp.bcf -- -t AN,AC",sep=""))
system(paste("mv /PROJECTS/ANONYMP/temp.bcf /PROJECTS/ANONYMP/chr15_5popSim_4B11_Ref.bcf",sep=""))
system(paste("tabix -f /PROJECTS/ANONYMP/chr15_5popSim_4B11_Ref.bcf",sep=""))

system(paste("/PROGS/EXTERN/shapeit/SHAPEIT5/shapeit5/xcftools/bin/xcftools view --input /PROJECTS/ANONYMP/chr15_5popSim_4B11_Ref.bcf -O sh -o /PROJECTS/ANONYMP/chr15_5popSim_Ref_4B11_xcf.bcf --r ",chr," --m 0",sep=""))

system(paste("bcftools convert -O b -o /PROJECTS/ANONYMP/chr15_5popSim_4B11_Cases.bcf /PROJECTS/ANONYMP/chr15_5popSim_4B11_Cases.vcf.gz",sep=""))
system(paste("tabix -f /PROJECTS/ANONYMP/chr15_5popSim_4B11_Cases.bcf",sep=""))

system(paste("bcftools +fill-tags /PROJECTS/ANONYMP/chr15_5popSim_4B11_Cases.bcf -Ob -o /PROJECTS/ANONYMP/temp.bcf -- -t AN,AC",sep=""))
system(paste("mv /PROJECTS/ANONYMP/temp.bcf /PROJECTS/ANONYMP/chr15_5popSim_4B11_Cases.bcf",sep=""))
system(paste("tabix -f /PROJECTS/ANONYMP/chr15_5popSim_4B11_Cases.bcf",sep=""))

system(paste("/PROGS/EXTERN/shapeit/SHAPEIT5/shapeit5/xcftools/bin/xcftools view --input /PROJECTS/ANONYMP/chr15_5popSim_4B11_Cases.bcf -O sh -o /PROJECTS/ANONYMP/chr15_5popSim_4B11_Cases_xcf.bcf --r ",chr," --m 0",sep=""))

ref1<-paste0("/PROJECTS/ANONYMP/chr15_5popSim_Ref_4B11_xcf")
target<-paste("/PROJECTS/ANONYMP/chr15_5popSim_4B11_Cases",sep="")

system(paste("/PROGS/EXTERN/IMPUTE5/impute5_v1.2.0/imp5Chunker_v1.2.0_static --h ",ref1,".bcf --g ",target,"_xcf.bcf --r ",chr," --o /PROJECTS/ANONYMP/coordinates.I51.2_",chr,".txt",sep=""))




####
###
##
####
###
##
corrdinates<-read.table(paste("/PROJECTS/ANONYMP/coordinates.I51.2_",chr,".txt",sep=""),header=F,as.is=T)

for (j in 1:nrow(corrdinates)){

system(paste("/PROGS/EXTERN/IMPUTE5/impute5_v1.2.0/impute5_static --h ",ref1,".bcf --g ",target,"_xcf.bcf --r ",corrdinates[j,4]," --m /PUBLIC_DATA/GeneticMaps/Bherer/Refined_EUR_genetic_map_b37/sexavg_chr15_s2.txt --surfbat --buffer-region ",corrdinates[j,3]," --o ",target,".I51.2_chunk",corrdinates[j,1],".vcf.gz --l ",target,".I51.2_chunk",corrdinates[j,1],".log",sep=""))

if(j==1){
write(paste(target,".I51.2_chunk",corrdinates[j,1],".vcf.gz",sep=""),paste("/PROJECTS/ANONYMP/ligate.I51.2_chr",chr,".txt",sep=""))
} else {
write(paste(target,".I51.2_chunk",corrdinates[j,1],".vcf.gz",sep=""),paste("/PROJECTS/ANONYMP/ligate.I51.2_chr",chr,".txt",sep=""),append=T)

}
}
####
###
##

system(paste("bcftools concat -n -f/PROJECTS/ANONYMP/ligate.I51.2_chr",chr,".txt -Oz -o ",target,".I51.2_chunkAll.vcf.gz",sep=""))

system(paste("tabix -f ",target,".I51.2_chunkAll.vcf.gz",sep=""))

for (j in 1:nrow(corrdinates)){
system(paste("rm ",target,".I51.2_chunk",corrdinates[j,1],".vcf.*",sep=""))
system(paste("rm ",target,".I51.2_chunk",corrdinates[j,1],".log",sep=""))
} 

system(paste("rm /PROJECTS/ANONYMP/ligate.I51.2_chr",chr,".txt",sep=""))






ref<-paste("/PROJECTS/ANONYMP/chr15_5popSim_4B11_Ref.vcf.gz",sep="")

system(paste("Minimac3 --refHaps ",ref," --processReference --prefix /PROJECTS/ANONYMP/chr15_5popSim_4B11_Ref.MinimacFull --chr ",chr,sep=""))

system(paste("/PROGS/EXTERN/Minimac3/Minimac4/release-build/minimac4 --update-m3vcf /PROJECTS/ANONYMP/chr15_5popSim_4B11_Ref.MinimacFull.m3vcf.gz > /PROJECTS/ANONYMP/chr15_5popSim_4B11_Ref.MinimacFull.msav --threads 16",sep=""))

system(paste("/PROGS/EXTERN/Minimac3/Minimac4/release-build/minimac4 -O vcf.gz -o /PROJECTS/ANONYMP/chr15_5popSim_4B11_Cases.mini4.vcf.gz -f GT,DS,HDS,SD,GP /PROJECTS/ANONYMP/chr15_5popSim_4B11_Ref.MinimacFull.msav /PROJECTS/ANONYMP/chr15_5popSim_4B11_Cases.vcf.gz --threads 16",sep=""))
system(paste("tabix -f /PROJECTS/ANONYMP/chr15_5popSim_4B11_Cases.mini4.vcf.gz",sep=""))




##########bash commands to extract dosages for making plot late on
##########


tmpdir=./


label=I5_anonymp

invcf=chr15_5popSim_4B11_Cases.I51.2_chunkAll.vcf.gz


zcat ${invcf} | grep -v "^#" | cut -f 1-5,8 --output-delimiter=" " | sed --regexp-extended 's/IMP;//g;s/AF=//g;s/INFO=//g;s/\;/ /g;s/AC=//g;s/AN=//g' | cut -d" " -f1-9 > ${tmpdir}SAP.${label}info.temp

k1=`zcat ${invcf} | grep -m1 "^#CHROM" | wc | awk '{print $2}'`
k1=$((5*($k1-9)))
seq 3 5 $k1 > ${tmpdir}cols.${label}.txt

zcat ${invcf} | grep -v "^#" | cut -f 10- --output-delimiter=" " | tr "[ ]" "[:]" | cut -d":" -f$(echo `sort -n ${tmpdir}cols.${label}.txt` | tr "[ ]" "[,]") | tr "[,:]" "[  ]" > ${tmpdir}SAP.${label}.temp

rm ${tmpdir}cols.${label}.txt




tmpdir=./


label=M4_anonymp

invcf=/PROJECTS/ANONYMP/chr15_5popSim_4B11_Cases.mini4.vcf.gz

zcat ${invcf} | grep -v "^#" | cut -f 1-5,8 --output-delimiter=" " | sed --regexp-extended 's/IMPUTED;//g;s/MAF=//g;s/R2=//g;s/\;/ /g;s/AF=//g;s/AVG_CS=//g' | cut -d" " -f1-9 > ${tmpdir}SAP.${label}info.temp

k1=`zcat ${invcf} | grep -m1 "^#CHROM" | wc | awk '{print $2}'`
k1=$((5*($k1-9)))
seq 2 5 $k1 > ${tmpdir}cols.${label}.txt

zcat ${invcf} | grep -v "^#" | cut -f 10- --output-delimiter=" " | tr "[ ]" "[:]" | cut -d":" -f$(echo `sort -n ${tmpdir}cols.${label}.txt` | tr "[ ]" "[,]") | tr "[,:]" "[  ]" > ${tmpdir}SAP.${label}.temp

rm ${tmpdir}cols.${label}.txt



###fin