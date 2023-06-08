
runSLCMA.adversity.20200406 <- function(adversity){
  
if(!adversity %in% c("abuse","faminstability","Fscore","mompsych","nbhqual","oneadult","parcruelty")){ 
  stop("Please specify an adversity measure from abuse, faminstability, Fscore, mompsych, nbhqual, oneadult, parcruelty")} 

cat("Loading packages \n")
suppressMessages( library(dplyr, lib.loc = "/data/js95/DunnLab/Yiwen/R_package/") )
suppressMessages( library(glmnet, lib.loc = "/data/js95/DunnLab/Yiwen/R_package/") )
suppressMessages( library(lars))
suppressMessages( library(intervals, lib.loc = "/data/js95/DunnLab/Yiwen/R_package/") )
suppressMessages( library(selectiveInference, lib.loc = "/data/js95/DunnLab/Yiwen/R_package/") )


# DNA methylation data and covariates prep -----------------------------------------------

cat("Loading DNA methylation data \n")
load("/data/js95/DunnLab/alussier/DNAm_samples_cells_15up_notwins.Rdata")
dim(betas.f15.win.notwins)
dim(samples.f15.notwins)
dim(cell.f15.notwins)

cat("Loading covariate data \n")
load("/data/js95/DunnLab/alussier/Adversity_DNAm_f15/covariates.f15.20200224.Rdata")
dim(covariate.data.f15)


# Adversity data prep -----------------------------------------------------

cat("Loading adversity data \n")
load("/data/js95/DunnLab/alussier/Adversity_DNAm_f15/adversity.data.f15.2020303.Rdata")
if(adversity %in% c("parcruelty", "Fscore")){
  print("Loading from KA.data")
  load("~/ALSPACnotwins.2018.08.27.Rdata")
  ka.data <- df
  dim(ka.data)
  ka.data <- ka.data[which(ka.data$cidB1471 %in% samples.f15.notwins$ALN),]
  ka.data <- ka.data[match(samples.f15.notwins$ALN, ka.data$cidB1471),]
  adversity.data.2 <- ka.data[,c(1, grep(adversity, colnames(ka.data)))]
  adversity.data.2 <- adversity.data.2[, -grep("18y", colnames(adversity.data.2))]
  rm(ka.data)
}

else{
  print("Loading from Beast")
  adversity.data.2 <- adversity.data[,c(1, grep(adversity, colnames(adversity.data)))] 
  }

rm(adversity.data)

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

#adversity.data.2$ever <- ifelse(adversity.data.2$accumulation==0, 0, 1)

save(adversity.data.2, file = paste("/data/js95/DunnLab/alussier/Adversity_DNAm_f15/", adversity,"/",adversity, "_exposure_data_", Sys.Date(),".Rdata", sep =""))

# SLCMA prep --------------------------------------------------------------
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

cat("Recoding to factors for the following variables: \n")
to.recode <- c(adversity.names, "ever","Female","WHITE","ed_momgest","mom_birthage","sustained.smoke")
cat(to.recode, "\n")
for(i in to.recode){
  lars.df[,i] <- as.factor(lars.df[,i])
}
#str(lars.df[,to.recode])

colnames(lars.df)[1] <- "ID"
sum(complete.cases(lars.df))

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


filename <- paste("/data/js95/DunnLab/alussier/Adversity_DNAm_f15/", adversity,"/", "LARS_", 
                  adversity, "_stage1_FWL_sI_noever_", Sys.Date(), ".Rdata",
                  sep =""
                   )
cat("Saving file as : \n", filename, "\n")

save(lars.res, file = filename)

cat("DONE \n")
}





