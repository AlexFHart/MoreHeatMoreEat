##					libraries						####


#library(Hmisc)
#library(car)
library(edgeR)
#library(tmod)
#library("ontologyIndex")
#library(qvalue)
#library(RRHO)
#library(lattice)
#help(readDGE)
#help(read.delim)
#library(vioplot)
#library(RColorBrewer)
#library (gplots)
#library(VennDiagram)


###					data						####
#############################################################################
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/PostDoc Folder/Scripts & Code/ThermAdapt/VectorAnalysis")
setwd("~/Documents/Alex/github/ThermAdapt/VectorAnalysis")

files<-read.table("Files_Info.tsv", stringsAsFactors=T, header=T)
table(files[,6])	#treatment
table(files[,5])	#experimental day
table(files[,4])	#time since mating: short, medium, long

counts<- readDGE(files[,1], path="~/Library/Mobile Documents/com~apple~CloudDocs/Documents/PostDoc Folder/HS_MR/Counts (new)", header=FALSE, 
	group=files[,6],labels=files[,2])
counts$samples

counts<- readDGE(files[,1], path="~/Documents/Alex/Thermal Niche Evo data/Vector Analysis/new_counts", header=FALSE, 
                 group=files[,6],labels=files[,2])


noint<-rownames(counts) %in% c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")
noint

counts<-counts[!noint,]
head(counts$counts)
tail(counts$counts)


counts$samples$lib.size <- colSums(counts$counts)


cpms <- cpm(counts)


############ filter everything that is not at least 3 cpm in at least one group

keep<-rowSums(cpms>3) >=5
length(keep[keep==TRUE])
length(keep[keep==FALSE])

counts<-counts[keep,]
#counts

counts$samples$lib.size <- colSums(counts$counts)

summary(counts$samples$lib.size/1000000)
counts$samples$lib.size



TestedGenes<-row.names(counts$counts)
length(TestedGenes)
#write.table(TestedGenes, "All_Genes_HS_Response.txt", quote=F, row.names=F,
#	col.names=F, sep="\t")

counts <- calcNormFactors(counts)

plotMDS(counts, 
	labels=counts$samples$group,
        col=as.numeric(factor(files$Time)),
	dim.plot=c(1,2))


short <- which(files[,4]=="S")
medium <- which(files[,4]=="M")
long <- which(files[,4]=="L")

plotMDS(counts[,short], 
        labels=counts[,short]$samples$group,
        col=as.numeric(factor(files[short,]$Treatment)),
        dim.plot=c(1,2))

plotMDS(counts[,medium], 
        labels=counts[,medium]$samples$group,
        col=as.numeric(factor(files[medium,]$Treatment)),
        dim.plot=c(1,2))

plotMDS(counts[,long], 
        labels=counts[,long]$samples$group,
        col=as.numeric(factor(files[long,]$Treatment)),
        dim.plot=c(1,2))

##					DEG analysis general				####

counts$samples
head(files)


treat <-files[,3]
time <-files[,4]
day <- files[,5]



treat<-relevel(treat, ref="FST")
time<-relevel(time, ref="S")
day<-as.factor(day)

var_matrix <- data.frame (treat=treat, day = day, time=time)


cbind(var_matrix, as.vector(counts$samples$group))



#-----------------------------------------------------------------------#

##												#
##												#
##					C complete data set				#
##												#
##												#

#-----------------------------------------------------------------------#


files$Group

design_C <-model.matrix(~ 0 + files$Group)
rownames(design_C) <- colnames(counts)
design_C

head(design_C)

v_C<-voom(counts,design_C, normalize.method="none", plot=TRUE)
fit_C <- lmFit(v_C, design_C)

head(v_C$E)


cont.matrix_C <- cbind(FEM_short=c(0,0,1,0,0,-1,0,0,0,0,0,0,0,0,0),
				HES_short=c(0,0,0,0,0,-1,0,0,1,0,0,0,0,0,0),
				IRR_short=c(0,0,0,0,0,-1,0,0,0,0,0,1,0,0,0),
				TRD_short=c(0,0,0,0,0,-1,0,0,0,0,0,0,0,0,1),
				FEM_medium=c(0,1,0,0,-1,0,0,0,0,0,0,0,0,0,0),
				HES_medium=c(0,0,0,0,-1,0,0,1,0,0,0,0,0,0,0),
				IRR_medium=c(0,0,0,0,-1,0,0,0,0,0,1,0,0,0,0),
				TRD_medium=c(0,0,0,0,-1,0,0,0,0,0,0,0,0,1,0),
				FEM_long=c(1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0),
				HES_long=c(0,0,0,-1,0,0,1,0,0,0,0,0,0,0,0),
				IRR_long=c(0,0,0,-1,0,0,0,0,0,1,0,0,0,0,0),
				TRD_long=c(0,0,0,-1,0,0,0,0,0,0,0,0,1,0,0),
				FEM = c(1,1,1,-1,-1,-1,0,0,0,0,0,0,0,0,0)/3,
				HES=c(0,0,0,-1,-1,-1,1,1,1,0,0,0,0,0,0)/3,
				IRR=c(0,0,0,-1,-1,-1,0,0,0,1,1,1,0,0,0)/3,
				TRD=c(0,0,0,-1,-1,-1,0,0,0,0,0,0,1,1,1)/3,
				medium_short=c(0,1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1)/5,
				long_medium=c(1,-1,0,1,-1,0,1,-1,0,1,-1,0,1,-1,0)/5,
				long_short=c(1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,0,-1)/5,
				Time_fst=c(0,0,0,1,0,-1,0,0,0,0,0,0,0,0,0),
				Time_trd=c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,-1))

fit2_C<- contrasts.fit(fit_C, cont.matrix_C)
fit3_C <- eBayes(fit2_C)

s_all <- summary(decideTests(fit3_C, method="separate", p.value=0.05))
s_all


#---------------------------------------------------------------------------

FEM_Top <-topTable(fit3_C, adjust="BH", coef="FEM",n=Inf, p=0.05, sort.by="logFC")
FEM_Complete<-topTable(fit3_C, adjust="BH", coef="FEM",n=Inf, p=1, sort.by="logFC")
hist(FEM_Complete$P.Value)

str(FEM_Top)
head(FEM_Top)

write.table(FEM_Top, "new_HS_vector_d2-3.txt", quote=F, sep="\t")


