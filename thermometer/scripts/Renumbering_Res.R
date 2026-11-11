###############################
#### renumber PDB files #######
###############################

library(bio3d)



pdb_format<-function(x, name) {
  
  matri<-x[,9:11]
  
  coord<-as.vector(apply(matri, 1, function(x) list=c(x[1], x[2], x[3])))
  
  write.pdb(xyz = coord, type=as.character(x$type), resno=as.character(x$resno), resid=as.character(x$resid), eleno=as.character(x$eleno), elety=as.character(x$elety), chain=as.character(x$chain), insert=as.character(x$insert),o = as.character(x$o), b = as.character(x$b), file=name)
  
}



##### read input


args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 2){
  
  print("Error!")
  print("Invalid number of arguments!")
  print("Run:")
  print("Rscript Renumbering_res.R dict_tmp.csv pathtoNewPdbFile")
  quit(status = 1)
}


dict_name <- args[1]
pdb_name_new <- args[2]



#### read the dictionary
dict <- read.csv(dict_name)



### resd the pdb structure file
struct <- read.pdb(pdb_name_new)
struct <- struct$atom
struct_2 <- struct

for(i in 1:nrow(dict)){
  
  struct_2[struct$resno == dict[i,2],]$chain <- strsplit(as.character(dict[i,1])," ")[[1]][1]
  struct_2[struct$resno == dict[i,2],]$resno <- strsplit(as.character(dict[i,1])," ")[[1]][2]
  if(length(strsplit(as.character(dict[i,1])," ")[[1]]) == 3){
    struct_2[struct$resno == dict[i,2],]$insert <- strsplit(as.character(dict[i,1])," ")[[1]][3]
  }
}


pdb_format(struct_2, pdb_name_new)