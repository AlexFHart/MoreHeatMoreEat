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
#library(statmod)

##						data						####

setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Documents/PostDoc Folder/Scripts & Code/ThermAdapt/VectorAnalysis")

files<-read.table("Files_Counts")
head(files)
files
lome <- which(files$V3=="L"|files$V3=="V")
files_lome <- files[lome,]

counts<- readDGE(files_lome[,5], path="~/Library/Mobile Documents/com~apple~CloudDocs/Documents/PostDoc Folder/HS_MR/Counts (new)", header=FALSE, 
	group = files_lome[,2], labels=files_lome[,2])
counts$samples



noint<-rownames(counts) %in% c("__no_feature", "__ambiguous", "__too_low_aQual", "__not_aligned", "__alignment_not_unique")
noint

counts<-counts[!noint,]
head(counts$counts)
tail(counts$counts)


counts$samples$lib.size <- colSums(counts$counts)


cpms <- cpm(counts)
head(cpms)
counts
str(cpms)
str(counts)
names(counts)
str(counts$samples)
counts$samples$group


############ filter everything that is not at least 1 cpm

keep<-rowSums(cpms>1) >=2
keep
length(keep[keep==TRUE])
length(keep[keep==FALSE])

counts<-counts[keep,]
counts

counts$samples$lib.size <- colSums(counts$counts)

summary(counts$samples$lib.size/1000000)
counts$samples$lib.size



TestedGenes<-row.names(counts$counts)
length(TestedGenes)
#write.table(TestedGenes, "All_Genes_Mating_Response.txt", quote=F,
#	sep="\t", row.names=F, col.names=F)

counts <- calcNormFactors(counts)

plotMDS(counts, 
	labels=counts$samples$group,
	col=c("red","blue", "blue","green", "black", "black", "red", "blue", "black", "green", "green"),
	dim.plot=c(1,2))
legend(0, 0.6, legend=c("poly","male_lim", "mono", "vrigin"), fill=c("blue", "red", "green", "black"))



##					DEG analysis general				####

counts$samples


treat <-factor(c("male_lim", "poly", "poly", "mono", "virgin", "virgin", "male_lim", "poly", "virgin",
                "mono", "mono"))


treat<-relevel(treat, ref="virgin")


var_matrix <- data.frame (treat=treat)


cbind(var_matrix, as.vector(counts$samples$group))

summary(var_matrix)


design_l <-model.matrix(~  0+treat)
rownames(design_l) <- colnames(counts)
design_l


plotMDS(counts, dim=c(1,2))


v_l<-voom(counts,design_l, normalize.method="none", plot=TRUE)



fit_l <- lmFit(v_l, design_l)


cont.matrix_l <- cbind(mating=c(-1, 1/3, 1/3, 1/3),	#overall mating response
				poly=c(-1,0,0,1),			#mating response to polygamous males
                        male_limited = c(-1,1,0,0),	#mating response to male-limited males
                        mono = c(-1,0,1,0),		#mating response to monogamous males
				mono_poly=c(0,0,1,-1),
				mal_poly=c(0,1,0,-1),
				mal_mono=c(0,1,-1,0))

fit2_l <- contrasts.fit(fit_l, cont.matrix_l)
fit3_l <- eBayes(fit2_l)

s_all <- summary(decideTests(fit3_l, method="separate", p.value=0.05))
s_all

s_test <- decideTests(fit3_l, method="separate", p.value=0.05)
str(s_test)
head(s_test)

#### genes that respond sign. to mating in all three settings (poly, mono & male lim males) ###

selectgenes <- s_test[,2]!=0 & s_test[,3]!=0 & s_test[,4]!=0
sum(selectgenes)

common_genes <- rownames(s_test[selectgenes,])

vennDiagram(vennCounts(s_test[,2:4], include="both"))


### effect sizes etc. for overall mating response

mat_all <- topTable(fit3_l, adjust="BH", coef=1,n=Inf, p=1, sort.by="logFC")
head(mat_all)

### filter for 1200 genes that respond sign. in all three settings 

top1200 <- mat_all[common_genes,]
head(top1200)

write.table(top1200, "new_MR_vector.txt", quote=F, sep="\t")
