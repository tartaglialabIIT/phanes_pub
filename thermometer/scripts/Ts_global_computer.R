########### This script computes the Ts score of a given protein structure #########

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 3){
  
  print("Error!")
  print("Invalid number of arguments!")
  print("Run:")
  print("Rscript Ts_computer.R EneFile Nthreads name_pdb")
  quit(status = 1)
}


### THERMO_HOME: root of the thermometer tool (scripts/, files/, bin/).
### Set so this script can run with cwd = a per-job work dir.
### Defaults to "." -> the original in-tool-dir invocation still works.
THERMO_HOME <- Sys.getenv("THERMO_HOME", unset = ".")

### Libraries

library(bio3d)
library(igraph)
#require(reshape2)

#require(stats)
#require(ape)
#require(plyr)
#require(stringr)
#require(foreach)
#require(doParallel)




EneFile <- args[1]
Nthread <- args[2]
name_pdb <- args[3]


    
#registerDoParallel(Nthread)
#getDoParWorkers()


Energy_Lim <- 150


### working directory definitions
# cwd is the per-job work dir; tool assets come from THERMO_HOME

#### Importing functions...
source(file.path(THERMO_HOME, "scripts", "Functions.R"))




### Importing Energy distributions...
df_density_range_lj <- read.csv(file.path(THERMO_HOME, "files", "df_density_range_lj.csv"))
df_density_range_C <- read.csv(file.path(THERMO_HOME, "files", "df_density_range_Coulomb.csv"))
df_density_range_E <- read.csv(file.path(THERMO_HOME, "files", "df_density_range_E_tot.csv"))



### Importing Energy names...
FileList_ene_csv <- c(EneFile)




################################################################################################################################################################
################################################################################################################################################################
################## 
################## Randomizing networks  with all Energy Terms (LJ , Coulomb and total)  in order to study  Strength parameter
################## 
################################################################################################################################################################
################################################################################################################################################################



#foreach(i=1:length(FileList_ene_csv)) %dopar% {
for(i in 1:length(FileList_ene_csv)){  
  
  vet_strength_TotE_Neg <- c()
  vet_strength_TotE_Pos <- c()
  vet_strength_TotE_Neg_Pos <- c()
  
  vet_strength_LJ_Neg <- c()
  vet_strength_LJ_Pos <- c()
  vet_strength_LJ_Neg_Pos <- c()
  
  vet_strength_C_Neg <- c()
  vet_strength_C_Pos <- c()
  vet_strength_C_Neg_Pos <- c()
  
  df_Strength <- data.frame()
  
  #setwd(dir_input)
  
  name <- as.character(FileList_ene_csv[i])
  print(name)
  df_ene <- read.csv(name)
  df_ene <- df_ene[df_ene$vet_SumE<150,c("vet_n","vet_contacts","vet_dist","vet_SumE_c","vet_SumE_lj","vet_SumE")]
  
  
  ######################## Tot Energy  ##################################################### 
  
  ### TRUE STRENGTH... 
  
  ave_strength_true_tot <- strength_fn(df_ene$vet_SumE, df_ene)
  ave_strength_true_C <- strength_fn(df_ene$vet_SumE_c, df_ene)
  ave_strength_true_lj <- strength_fn(df_ene$vet_SumE_lj, df_ene)
  
  
  ### RANDOM STRENGTH...
  print("Randomization of total energy...")
  ene_rand_aus_tot <- t(apply(df_ene[3],1,rand_ene,df_density_range_E))
  
  print("Randomization of coulomb energy...")
  ene_rand_aus_C <- t(apply(df_ene[3],1,rand_ene,df_density_range_C))
  
  print("Randomization of lj energy...")
  ene_rand_aus_lj <- t(apply(df_ene[3],1,rand_ene,df_density_range_lj))
  
  
  print("Computing strength on total energy...")
  strength_rand_tot <- apply(ene_rand_aus_tot, 2, strength_fn, df_ene)
  
  print("Computing strength on coulomb energy...")
  strength_rand_C <- apply(ene_rand_aus_C, 2, strength_fn, df_ene)
  
  print("Computing strength on lj energy...")
  strength_rand_lj <- apply(ene_rand_aus_lj, 2, strength_fn, df_ene)
  
  
  ### Strength vectors...
  vet_strength_tot <- cbind(ave_strength_true_tot, strength_rand_tot)
  vet_strength_C <- cbind(ave_strength_true_C, strength_rand_C)
  vet_strength_lj <- cbind(ave_strength_true_lj, strength_rand_lj)
  
  ########################################################################################
  
  
  ### Saving Strength values for each protein
  df_Strength <- t(rbind(vet_strength_tot,vet_strength_C, vet_strength_lj))
  
  colnames(df_Strength) <- c("S_tot_PN","S_tot_N", "S_tot_P", "S_C_PN","S_C_N", "S_C_P", "S_lj_PN","S_lj_N", "S_lj_P")
  
  #n <- unlist(strsplit(name, "[.]"))[1]
  file_name <- paste(name_pdb, "Strength_Global.csv", sep="_")
  
  
  write.csv(df_Strength, file = file_name)
  
}


###############  Computing Ts descriptors...




#Network_Parameter <- "Strength_Global.csv"
#FileList_rand_csv <- list.files(dir_output, Network_Parameter)
#print(FileList_rand_csv)

FileList_rand_csv <- c(file_name)

###### ROC and PERFORMANCE ANALYSIS #####
df_info <- data.frame()
prot_vet <- c()

for (i in 1:length(FileList_rand_csv)) {
  
  name <- as.character(FileList_rand_csv[i])
  df_rand <- read.csv(name)
  
  # Use the id the caller passed, not strsplit(filename,"_")[1] -- the old
  # form turned "AF-O00168-F1-model_v6_clean" into "AF-O00168-F1-model",
  # which no longer matched the keys used by the other feature tables.
  prot <-  name_pdb
  
  d <- dim(df_rand)[2]
  len <- dim(df_rand)[1]
  vet_ave_S <- c()
  vet_sd_S <- c()
  
  for (j in 2:d) {
    average_S <- mean(as.numeric(as.vector(df_rand[2:len,j])))
    sd_S <- sd(as.numeric(as.vector(df_rand[2:len,j])))
    vet_ave_S <- c(vet_ave_S, average_S)
    vet_sd_S <- c(vet_sd_S, sd_S)
  }
  
  vet_ave_true <- df_rand[1,2:d]
  
  prot_vet <- c(prot_vet, prot)
  df_info <- rbind(df_info, c(as.numeric(as.vector(vet_ave_true)),vet_ave_S, vet_sd_S))
  
}

df_res <- as.data.frame(cbind(prot_vet, df_info))
dim(df_res)

colnames(df_res) <- c("PDB", "S_tot_PN_Real","S_tot_N_Real", "S_tot_P_Real", 
                      "S_C_PN_Real","S_C_N_Real", "S_C_P_Real", "S_lj_PN_Real","S_lj_N_Real", 
                      "S_lj_P_Real","S_tot_PN_Mean","S_tot_N_Mean", "S_tot_P_Mean", "S_C_PN_Mean",
                      "S_C_N_Mean", "S_C_P_Mean", "S_lj_PN_Mean","S_lj_N_Mean", "S_lj_P_Mean", 
                      "S_tot_PN_sd","S_tot_N_sd", "S_tot_P_sd", "S_C_PN_sd",
                      "S_C_N_sd", "S_C_P_sd", "S_lj_PN_sd","S_lj_N_sd", "S_lj_P_sd")


write.csv(df_res, file = paste0(name_pdb, "_", "RealRandom_Strength.csv"))


Ts_score_tot_PN <- (df_res$S_tot_PN_Real - (df_res$S_tot_PN_Mean - df_res$S_tot_PN_sd))
Ts_score_tot_N <- (df_res$S_tot_N_Real - (df_res$S_tot_N_Mean - df_res$S_tot_N_sd))
Ts_score_tot_P <- (df_res$S_tot_P_Real - (df_res$S_tot_P_Mean - df_res$S_tot_P_sd))

Ts_score_C_PN <- (df_res$S_C_PN_Real - (df_res$S_C_PN_Mean - df_res$S_C_PN_sd))
Ts_score_C_N <- (df_res$S_C_N_Real - (df_res$S_C_N_Mean - df_res$S_C_N_sd))
Ts_score_C_P <- (df_res$S_C_P_Real - (df_res$S_C_P_Mean - df_res$S_C_P_sd))

Ts_score_lj_PN <- (df_res$S_lj_PN_Real - (df_res$S_lj_PN_Mean - df_res$S_lj_PN_sd))
Ts_score_lj_N <- (df_res$S_lj_N_Real - (df_res$S_lj_N_Mean - df_res$S_lj_N_sd))
Ts_score_lj_P <- (df_res$S_lj_P_Real - (df_res$S_lj_P_Mean - df_res$S_lj_P_sd))


df_Ts <- as.data.frame(cbind(prot_vet, Ts_score_tot_PN, Ts_score_tot_N, Ts_score_tot_P,  Ts_score_C_PN, Ts_score_C_N, Ts_score_C_P, Ts_score_lj_PN, Ts_score_lj_N, Ts_score_lj_P))

colnames(df_Ts) <- c("PDB", "Ts_tot_PN", "Ts_tot_N", "Ts_tot_P",  "Ts_C_PN", "Ts_C_N", "Ts_C_P", "Ts_lj_PN", "Ts_lj_N", "Ts_lj_P")


write.csv(df_Ts, file = paste0(name_pdb,"_","Ts_global_score.csv"))


