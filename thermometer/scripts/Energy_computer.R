########### This script computes the energy a given protein structure #########

args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 3){
  
  print("Error!")
  print("Invalid number of arguments!")
  print("Run:")
  print("Rscript Energy_computer.R pdbfile name_ene Nthreads")
  quit(status = 1)
}




### THERMO_HOME: root of the thermometer tool (scripts/, files/, bin/).
### Set so this script can run with cwd = a per-job work dir.
### Defaults to "." -> the original in-tool-dir invocation still works.
THERMO_HOME <- Sys.getenv("THERMO_HOME", unset = ".")

### Libraries
library(igraph)
library(reshape2)
library(bio3d)
library(svMisc)


#require(stats)
#require(ape) 
#require(plyr) 
#require(pvclust)
#require(foreach)
#require(doParallel)




pdbfile <- args[1]
name_ene <- args[2]
Nthread <- args[3]



#registerDoParallel(Nthread)
#getDoParWorkers()

### Working directories
# cwd is the per-job work dir; tool assets come from THERMO_HOME

# Loading required functions... 
source(file.path(THERMO_HOME, "scripts", "Functions.R"))



LJcoeff <- read.csv(file.path(THERMO_HOME, "files", "ljcoeff.txt"), header = T, sep = ",")
atomMap <- read.csv(file.path(THERMO_HOME, "files", "Cariche.txt"), header = T, sep = ",")

### Upload a list of PDB in the input directory
ListPDB <- c(pdbfile)
print(paste0("Length of list of pdb -->", length(ListPDB)))
### You can run one of them or run them all (Current choise <- run them all)

#foreach(k=1:length(ListPDB)) %dopar% {
for(k in 1:length(ListPDB)){
    
  # Initialize DF
  PDB_Energy <- data.frame()
  
  # Upload PDB and create elechain column (ID of single residue)
  pdb <- read.pdb(as.character(ListPDB[k]))
  pdb$atom[is.na(pdb$atom$insert),8] <- ""
  
  
  # Run Energy function
  PDB_Energy <- energies(pdb, LJcoeff, atomMap)
  
  # Save csv file
  #pdb_name <- sub("*_min.pdb", "", ListPDB[k])
  write.csv(PDB_Energy, file = paste0(name_ene, '.csv'), quote = FALSE, eol = "\n", na = "NA", row.names = TRUE, fileEncoding = "")
  
}
