########### This script computes the Ts score at single residue of a given protein structure #########

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 3){
  
  print("Error!")
  print("Invalid number of arguments!")
  print("Run:")
  print("Rscript Ts_local_computer.R ene_pdb  Nthreads pdb_name")
  quit(status = 1)
}



### THERMO_HOME: root of the thermometer tool (scripts/, files/, bin/).
### Set so this script can run with cwd = a per-job work dir.
### Defaults to "." -> the original in-tool-dir invocation still works.
THERMO_HOME <- Sys.getenv("THERMO_HOME", unset = ".")

### Libraries
library(bio3d)
library(igraph)

#require(foreach)
#require(doParallel)



ene_pdb <- args[1]
Nthread <- args[2]
pdb_name <- args[3]



#registerDoParallel(Nthread)
#getDoParWorkers()


Energy_Lim <- 150

# cwd is the per-job work dir; tool assets come from THERMO_HOME


#### Importing functions...
source(file.path(THERMO_HOME, "scripts", "Functions.R"))


### Importing Energy distributions...
df_density_range_lj <- read.csv(file.path(THERMO_HOME, "files", "df_density_range_lj.csv"))
df_density_range_C <- read.csv(file.path(THERMO_HOME, "files", "df_density_range_Coulomb.csv"))
df_density_range_E <- read.csv(file.path(THERMO_HOME, "files", "df_density_range_E_tot.csv"))




### Importing Energy names...
FileList_ene_csv <- c(ene_pdb) 

ResMatrix <- c()

################################################################################################################################################################
################################################################################################################################################################
################## 
################## Randomizing networks  with all Energy Terms (LJ , Coulomb and total)  in order to study  Strength parameter
################## 
################################################################################################################################################################
################################################################################################################################################################

#Total_graph_par <- data.frame()

#foreach(i=1:length(FileList_ene_csv)) %dopar% {
for(i in 1:length(FileList_ene_csv)){
  
  
  
  
  name <- as.character(FileList_ene_csv[i])
  #print(name)
  df_ene <- read.csv(name)
  df_ene <- df_ene[df_ene$vet_SumE<Energy_Lim,c("vet_n","vet_contacts","vet_dist","vet_SumE_c","vet_SumE_lj","vet_SumE")]
  
  #setwd(dir_function)
  #source("random_function.R")
  
  
  ######################## Tot Energy  ##################################################### 
  
  
  
  ### RANDOM STRENGTH...
  print("Randomization of total energy...")
  ene_rand_aus_tot <- t(apply(df_ene[3],1,rand_ene,df_density_range_E))
 
  print("Randomization of coulomb energy...")
  ene_rand_aus_C <- t(apply(df_ene[3],1,rand_ene,df_density_range_C))
  
  print("Randomization of lj energy...")
  ene_rand_aus_lj <- t(apply(df_ene[3],1,rand_ene,df_density_range_lj))
  
  
  
  print("Computing strength on total energy...")
  strength_true_tot <- strength_Spin_fn(df_ene$vet_SumE, df_ene)

  print("Computing strength on couplomb energy...")
  strength_true_C <- strength_Spin_fn(df_ene$vet_SumE_c, df_ene)
  
  print("Computing strength on lj energy...")
  strength_true_lj <- strength_Spin_fn(df_ene$vet_SumE_lj, df_ene)
  
  print("E")
  strength_rand_tot <- apply(ene_rand_aus_tot, 2, strength_Spin_fn, df_ene)
  print("C")
  strength_rand_C <- apply(ene_rand_aus_C, 2, strength_Spin_fn, df_ene)
  print("lj")
  strength_rand_lj <- apply(ene_rand_aus_lj, 2, strength_Spin_fn, df_ene)
  
  ave_E_tot <- apply(strength_rand_tot, 1,mean)
  dev_E_tot <- apply(strength_rand_tot, 1,sd)
  
  ave_E_C <- apply(strength_rand_C, 1,mean)
  dev_E_C <- apply(strength_rand_C, 1,sd)
  
  ave_E_lj <- apply(strength_rand_lj, 1,mean)
  dev_E_lj <- apply(strength_rand_lj, 1,sd)
  
  Ts_tot <- strength_true_tot - (ave_E_tot - dev_E_tot) 
  Ts_C <-  strength_true_C - (ave_E_C - dev_E_C) 
  Ts_lj <- strength_true_lj - (ave_E_lj - dev_E_lj) 
  
  
  res <- as.data.frame(cbind(strength_true_tot, strength_true_C, strength_true_lj, ave_E_tot, dev_E_tot, ave_E_C, dev_E_C, ave_E_lj, dev_E_lj, Ts_tot, Ts_C, Ts_lj))
  
  colnames(res) <- c("S_true_tot", "S_true_C", "S_true_lj", "S_ave_tot", "S_dev_tot", "S_ave_C", "S_dev_C", "S_ave_lj", "S_dev_lj", "Ts_residue_tot", "Ts_residue_C", "Ts_residue_lj")
  filename <- paste(pdb_name, "Strength_SingleResidue.csv", sep = "_")
  write.csv(file = filename, x =res)
  
  
  
}



