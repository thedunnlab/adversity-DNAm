
setwd("/data/js95/DunnLab/alussier/Adversity_DNAm_f15")
#loading old data
load("/data/js95/DunnLab/alussier/Adversity_DNAm_f15/parcruelty/parcruelty_exposure_data_2020-05-26.Rdata", verbose=T)
old.parc <- adversity.data.2

load("parcrueltyRevised_2022-07-13.Rdata")
adversity = "parcruelty"
parcrueltyRevised$ID <- paste0(parcrueltyRevised$cidB1471, parcrueltyRevised$qlet)

suppressMessages( library(dplyr, lib.loc = "/data/js95/DunnLab/Yiwen/R_package/") )
suppressMessages( library(glmnet, lib.loc = "/data/js95/DunnLab/Yiwen/R_package/") )
suppressMessages( library(lars))
suppressMessages( library(intervals, lib.loc = "/data/js95/DunnLab/Yiwen/R_package/") )
suppressMessages( library(selectiveInference, lib.loc = "/data/js95/DunnLab/Yiwen/R_package/") )
library(selectiveInference)
cat("Loading DNA methylation data \n")
load("/data/js95/DunnLab/alussier/DNAm_samples_cells_15up_notwins.Rdata", verbose=T)
dim(betas.f15.win.notwins)
dim(samples.f15.notwins)
dim(cell.f15.notwins)

cat("Loading covariate data \n")
load("/data/js95/DunnLab/alussier/Adversity_DNAm_f15/covariates.f15.20200224.Rdata")
dim(covariate.data.f15)


parcrueltyRevised.2 <- parcrueltyRevised[which(parcrueltyRevised$ID %in% 
                                       paste0(samples.f15.notwins$ALN, samples.f15.notwins$QLET)),]
parcrueltyRevised.2 <- parcrueltyRevised.2[match(paste0(samples.f15.notwins$ALN, samples.f15.notwins$QLET), 
                                       parcrueltyRevised.2$ID),]
adversity.data.2 <- parcrueltyRevised.2[,c(1, grep(adversity, colnames(parcrueltyRevised.2)))]
adversity.data.2$cidB1471 <- samples.f15.notwins$ALN

adversity.names <- colnames(adversity.data.2)[-1]


cat("Making recency vector \n")
recency.vec <- gsub(paste(adversity, "_", sep =""), "", adversity.names)
recency.vec.2 <- recency.vec
recency.vec <- as.numeric(gsub("D","", 
                               gsub("y","", 
                                    gsub("m","", recency.vec))) )
recency.vec[grep("m", recency.vec.2)] <- recency.vec[grep("m", recency.vec.2)]/12

cat("Making accumulation and ever exposed variables \n")
adversity.data.2$accumulation <- rowSums(adversity.data.2[,adversity.names])
adversity.data.2$ever <- ifelse(rowSums(adversity.data.2[,adversity.names], na.rm =T) >0, 1, ifelse(rowSums(adversity.data.2[,adversity.names])==0, 0, NA))

save(adversity.data.2, file = paste("/data/js95/DunnLab/alussier/Adversity_DNAm_f15/parcruelty/",adversity, "_exposure_data_", Sys.Date(),".Rdata", sep =""))

summary(adversity.data.2)


cat("Preparing data for SLCMA \n")

sample.check <- identical(as.numeric(adversity.data.2$cidB1471), as.numeric(samples.f15.notwins$ALN))
if(sample.check ==F){stop("Samples in adversity data do not match DNA methylation data")}
sample.check <- identical(as.numeric(adversity.data.2$cidB1471), as.numeric(covariate.data.f15$cidB1471))
if(sample.check ==F){stop("Samples in adversity data do not match covariate samples")}

lars.df <- cbind(adversity.data.2, 
                 covariate.data.f15[,-1],
                 age = samples.f15.notwins$age,
                 betas.f15.win.notwins)
dim(lars.df)
sum(complete.cases(lars.df)) #661 - different from before

identical(adversity.data.2$cidB1471, old.parc$cidB1471)#TRUE
temp <- cbind(old.parc, 
              covariate.data.f15[,-1],
              age = samples.f15.notwins$age,
              betas.f15.win.notwins)
sum(complete.cases(temp)) #651

#trying to run with the same samples 
id.new <- lars.df$cidB1471[complete.cases(lars.df)]
id.old <- temp$cidB1471[complete.cases(temp)]
length(id.new[!id.new %in% id.old] ) #10
toremove <- id.new[!id.new %in% id.old] 
lars.df <- lars.df[!lars.df$cidB1471 %in% toremove, ]


cat("Recoding to factors for the following variables: \n")
to.recode <- c(adversity.names, "ever","Female","WHITE","ed_momgest","mom_birthage","sustained.smoke")
cat(to.recode, "\n")
for(i in to.recode){
  lars.df[,i] <- as.factor(lars.df[,i])
}

colnames(lars.df)[1] <- "ID"



# SLCMA -------------------------------------------------------------------
cat("Running SLCMA \n")
cat(date(), "\n")

covars <- c(colnames(covariate.data.f15[,-1]))
outcome.vec <- colnames(betas.f15.win.notwins)

source("/data/js95/DunnLab/alussier/LARS-noimpute-function-20190516.R")
lars.res <- select.LARS.complete(lars.df, 
                                 outcome.vec = outcome.vec, 
                                 adver = adversity,
                                 covars = covars, 
                                 hypos = c("accumulation","recency"),
                                 exposures <- 'default',
                                 recency.vec = recency.vec ,
                                 inf.method = "sI",
                                 FWL = T)  

filename <- paste("/data/js95/DunnLab/alussier/Adversity_DNAm_f15/parcruelty/", "LARS_", 
                  adversity, "_stage1_FWL_sI_noever_overlapOld_", Sys.Date(), ".Rdata",
                  sep ="")
cat("Saving file as : \n", filename, "\n")

save(lars.res, file = filename)

