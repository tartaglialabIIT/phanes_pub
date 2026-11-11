###############################
#### check PDB files ##########
###############################

library(bio3d)



pdb_format<-function(x, name) {
  
  matri<-x[,9:11]
  
  coord<-as.vector(apply(matri, 1, function(x) list=c(x[1], x[2], x[3])))
  
  write.pdb(xyz = coord, type=as.character(x$type), resno=as.character(x$resno), resid=as.character(x$resid), eleno=as.character(x$eleno), elety=as.character(x$elety), chain=as.character(x$chain), insert=as.character(x$insert),o = as.character(x$o), b = as.character(x$b), file=name)
  
}




### THERMO_HOME: root of the thermometer tool (scripts/, files/, bin/).
### Set so this script can run with cwd = a per-job work dir.
### Defaults to "." -> the original in-tool-dir invocation still works.
THERMO_HOME <- Sys.getenv("THERMO_HOME", unset = ".")

##### read input


args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 2){
  
  print("Error!")
  print("Invalid number of arguments!")
  print("Run:")
  print("Rscript Check_pdb.R pathtoPdbFile Max_number_residues")
  quit(status = 1)
}


pdb_name <- args[1]

# set the maximum allowed number of residues in the input PDB file 
thr_num_res <- as.numeric(args[2])


struct <- read.pdb(pdb_name)
struct <- struct$atom
struct <- struct[struct$type == "ATOM",]


if(length(struct$chain[is.na(struct$chain)]) == nrow(struct)){
  struct$chain <- "A"
  print("THERMOMETER output: No chain detected. All residues were labeled with the A chain identifier.", quote = F)
}


if(length(struct$chain[is.na(struct$chain)]) > 0 ){
  chain_number <- length(unique(struct$chain[! is.na(struct$chain)]))
  possible_chain <- LETTERS[! LETTERS %in% unique(struct$chain[! is.na(struct$chain)])]
  struct$chain[is.na(struct$chain)] <- possible_chain[1] 
  print(paste("THERMOMETER output: Some atoms are not labeled with a valid chain identifier. They are assigned to ", possible_chain[1], " chain.", sep=""),quote = F)
}




allowed_residues <- read.csv(file.path(THERMO_HOME, "files", "Cariche.txt"))
allowed_residues <- as.character(unique(allowed_residues$Resid))
allowed_residues <- c(allowed_residues, "HIS")

forbidden_residues <- struct[! struct$resid %in% allowed_residues,]
if(nrow(forbidden_residues) > 0 ){
  print("THERMOMETER output: Warning: Unkwown residue name(s) detected. Please check the PDB input file.",quote = F)
  print(paste("THERMOMETER output: Unknown residue names: ", unique(forbidden_residues$resid), sep=" "),quote = F)
  print("THERMOMETER output: ERROR: Please fix the input structure and resubmit.",quote = F)
  quit(status = 1)
}

struct <- struct[struct$resid %in% allowed_residues,]



struct_CA <- struct[struct$elety == "CA",]
stopifnot(nrow(struct_CA) <= thr_num_res)



chain_number <- unique(struct$chain[! is.na(struct$chain)])

for(i in 1:length(chain_number)){
  single_chain <- struct[struct$chain == chain_number[i],]
  if(length(table(diff(unique(single_chain$resno)) > 1)) > 1){
    print(paste("THERMOMETER output: Warning: Missing Residue(s) in chain ",chain_number[i],sep = ""),quote = F)
    print("THERMOMETER output: ERROR: Please fix the input structure and resubmit.",quote = F)
    quit(status = 1)
    
  }
}




struct[is.na(struct$insert),]$insert <- ""
residue_names <- unique(paste(struct$chain,struct$resno, struct$insert))
residue_chains <- unique(substring(residue_names,1,1))

dict_tot <- c()
for(p in 1:length(residue_chains)){
  tmp <- residue_names[substring(residue_names,1,1) == residue_chains[p]]
  dict <- cbind(tmp,seq(1,length(tmp),1))
  dict_tot <- rbind(dict_tot,dict)
  
}

write.csv(dict_tot, paste0(tools::file_path_sans_ext(pdb_name), "_dict_tmp.csv"), quote=F, row.names = F)


struct_2 <- struct
for(i in 1:nrow(dict_tot)){
  
  struct_2[paste(struct$chain,struct$resno, struct$insert) %in% dict_tot[i,1],]$resno <- dict_tot[i,2]
  struct_2[paste(struct$chain,struct$resno, struct$insert) %in% dict_tot[i,1],]$insert <- ""
}


pdb_format(struct_2, paste0(tools::file_path_sans_ext(pdb_name), "_clean.pdb"))

