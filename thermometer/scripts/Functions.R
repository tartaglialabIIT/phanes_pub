######## Number of random samble is set to:  1000

################## INIZIO IL PROCESSO DI RANDOMIZZAZIONE con casi fisici (energie estrapolate a determinati range di distanza)

rand_ene <- function(dist_link, df_density_range){
  
  Nsample = 1000
  row_density <- df_density_range[df_density_range$max_bin > dist_link,][1,]
  bin_size <- (row_density$max_ene - row_density$min_ene)/row_density$Num_Bin
  
  vet_mids <- c(0:(row_density$Num_Bin-1))*bin_size + row_density$min_ene + bin_size/2
  vet_prob <- row_density[,7:length(row_density)]*bin_size
  
  sample_row <- sample(vet_mids, vet_prob, size=Nsample, replace = T)
  #print("Qualche secondo di pazienza...")
  
  return(sample_row)
  
}


hubscore_fn <- function(ene_rand, df_ene_prot){
  
  df_aus <- as.data.frame(cbind(df_ene_prot[,c("vet_n", "vet_contacts")], ene_rand))
  
  EdgeList <- as.data.frame(df_aus)
  EdgeList[,1]=as.character(EdgeList[,1])
  EdgeList[,2]=as.character(EdgeList[,2])
  EdgeList <- as.matrix(EdgeList)
  g <- graph.edgelist(EdgeList[,1:2], directed=FALSE)
  E(g)$weight <- as.numeric(EdgeList[,3])
  
  Hub_Score <- hub_score(g, weights = E(g)$weight)
  Hub_Score <- Hub_Score$vector
  average_hubscore<- mean(Hub_Score)
  
}


strength_fn <- function(ene_rand, df_ene_prot){
  
  strenth_network_fn <- function(df_aus){
    
    EdgeList <- as.data.frame(df_aus)
    EdgeList[,1]=as.character(EdgeList[,1])
    EdgeList[,2]=as.character(EdgeList[,2])
    EdgeList <- as.matrix(EdgeList)
    g <- graph.edgelist(EdgeList[,1:2], directed=FALSE)
    E(g)$weight <- as.numeric(EdgeList[,3])
    
    Strength <- strength(g, weights = E(g)$weight)
    average_Strength <- mean(Strength)
    
    return(average_Strength)
  }
  
  
  
  df_tot <- as.data.frame(cbind(df_ene_prot[,c("vet_n", "vet_contacts")], ene_rand))
  colnames(df_tot) <- c("vet_n", "vet_contacts","ene_rand")
  
  average_strenth_tot <- strenth_network_fn(df_tot)
  average_strenth_neg <- strenth_network_fn(df_tot[df_tot$ene_rand <= 0, ])
  average_strenth_pos <- strenth_network_fn(df_tot[df_tot$ene_rand > 0, ])
  
  
  return(c(average_strenth_tot,average_strenth_neg,average_strenth_pos))
  
}




transitivity_fn <- function(ene_rand, df_ene_prot){
  
  transitivity_network_fn <- function(df_aus){
    
    EdgeList <- as.data.frame(df_aus)
    EdgeList[,1]=as.character(EdgeList[,1])
    EdgeList[,2]=as.character(EdgeList[,2])
    EdgeList <- as.matrix(EdgeList)
    g <- graph.edgelist(EdgeList[,1:2], directed=FALSE)
    E(g)$weight <- as.numeric(EdgeList[,3])
    
    Transitivity <- transitivity(g,  type="barrat", weights = E(g)$weight)
    average_Transitivity <- mean(Transitivity)
    
    return(average_Transitivity)
  }
  
  
  
  df_tot <- as.data.frame(cbind(df_ene_prot[,c("vet_n", "vet_contacts")], ene_rand))
  colnames(df_tot) <- c("vet_n", "vet_contacts","ene_rand")
  
  average_transitivity_tot <- transitivity_network_fn(df_tot)
  average_transitivity_neg <- transitivity_network_fn(df_tot[df_tot$ene_rand <= 0, ])
  average_transitivity_pos <- transitivity_network_fn(df_tot[df_tot$ene_rand > 0, ])
  
  
  return(c(average_transitivity_tot,average_transitivity_neg,average_transitivity_pos))
  
}


degree_fn <- function(ene_rand, df_ene_prot){
  
  degree_network_fn <- function(df_aus){
    
    EdgeList <- as.data.frame(df_aus)
    EdgeList[,1]=as.character(EdgeList[,1])
    EdgeList[,2]=as.character(EdgeList[,2])
    EdgeList <- as.matrix(EdgeList)
    g <- graph.edgelist(EdgeList[,1:2], directed=FALSE)
    E(g)$weight <- as.numeric(EdgeList[,3])
    
    Degree <-  degree(g, loops = FALSE, normalized = TRUE)
    average_Degree <- mean(Degree)
    
    return(average_Degree)
  }
  
  
  
  df_tot <- as.data.frame(cbind(df_ene_prot[,c("vet_n", "vet_contacts")], ene_rand))
  colnames(df_tot) <- c("vet_n", "vet_contacts","ene_rand")
  
  average_degree_tot <- degree_network_fn(df_tot)
  average_degree_neg <- degree_network_fn(df_tot[df_tot$ene_rand <= 0, ])
  average_degree_pos <- degree_network_fn(df_tot[df_tot$ene_rand > 0, ])
  
  
  return(c(average_degree_tot,average_degree_neg,average_degree_pos))
  
}



### working directory definitions
strength_Spin_fn <- function(ene_rand, df_ene_prot){
  
  strenth_network_fn <- function(df_aus){
    
    EdgeList <- as.data.frame(df_aus)
    EdgeList[,1]=as.character(EdgeList[,1])
    EdgeList[,2]=as.character(EdgeList[,2])
    EdgeList <- as.matrix(EdgeList)
    g <- graph.edgelist(EdgeList[,1:2], directed=FALSE)
    E(g)$weight <- as.numeric(EdgeList[,3])
    
    Strength <- strength(g, weights = E(g)$weight)
    #average_Strength <- mean(Strength)
    
    return(Strength)
  }
  
  
  
  df_tot <- as.data.frame(cbind(df_ene_prot[,c("vet_n", "vet_contacts")], ene_rand))
  colnames(df_tot) <- c("vet_n", "vet_contacts","ene_rand")
  
  strenth_tot <- strenth_network_fn(df_tot)
  #average_strenth_neg <- strenth_network_fn(df_tot[df_tot$ene_rand <= 0, ])
  #average_strenth_pos <- strenth_network_fn(df_tot[df_tot$ene_rand > 0, ])
  
  
  return(strenth_tot)#c(average_strenth_tot,average_strenth_neg,average_strenth_pos))
  
}



degree_Spin_fn <- function(ene_rand, df_ene_prot){
  
  degree_network_fn <- function(df_aus){
    
    EdgeList <- as.data.frame(df_aus)
    EdgeList[,1]=as.character(EdgeList[,1])
    EdgeList[,2]=as.character(EdgeList[,2])
    EdgeList <- as.matrix(EdgeList)
    g <- graph.edgelist(EdgeList[,1:2], directed=FALSE)
    E(g)$weight <- as.numeric(EdgeList[,3])
    
    Degree <-  degree(g, loops = FALSE, normalized = FALSE)
    #average_Degree <- mean(Degree)
    
    return(Degree)
  }
  
  
  
  df_tot <- as.data.frame(cbind(df_ene_prot[,c("vet_n", "vet_contacts")], ene_rand))
  colnames(df_tot) <- c("vet_n", "vet_contacts","ene_rand")
  
  average_degree_tot <- degree_network_fn(df_tot)
  #average_degree_neg <- degree_network_fn(df_tot[df_tot$ene_rand <= 0, ])
  #average_degree_pos <- degree_network_fn(df_tot[df_tot$ene_rand > 0, ])
  
  
  return(c(average_degree_tot))#,average_degree_neg,average_degree_pos))
  
}

self_strength_fn <- function(ene_rand, df_ene_prot){
  
  strenth_network_fn <- function(df_aus){
    
    EdgeList <- as.data.frame(df_aus)
    EdgeList[,1]=as.character(EdgeList[,1])
    EdgeList[,2]=as.character(EdgeList[,2])
    EdgeList <- as.matrix(EdgeList)
    g <- graph.edgelist(EdgeList[,1:2], directed=FALSE)
    E(g)$weight <- as.numeric(EdgeList[,3])
    
    Degree <- degree(g)
    
    Strength <- strength(g, weights = E(g)$weight)
    Strength <- Strength/Degree
    
    
    average_Strength <- mean(Strength)
    
    return(average_Strength)
  }
  
  
  df_tot <- as.data.frame(cbind(df_ene_prot[,c("vet_n", "vet_contacts")], ene_rand))
  colnames(df_tot) <- c("vet_n", "vet_contacts","ene_rand")
  
  average_strenth_tot <- strenth_network_fn(df_tot)
  average_strenth_neg <- strenth_network_fn(df_tot[df_tot$ene_rand <= 0, ])
  average_strenth_pos <- strenth_network_fn(df_tot[df_tot$ene_rand > 0, ])
  
  
  return(c(average_strenth_tot,average_strenth_neg,average_strenth_pos))
  
}






self_rand_ene <- function(df_ene){
  
  dr <- 0.05
  N_links <- dim(df_ene)[1]
  n <- 500
  
  df_res <- as.data.frame(df_ene[,4])
  
  for(j in 1:n){
    # print(j)
    df_res[,1+j] <- as.numeric(as.vector(df_ene[,4]))
  }  
  
  
  len <- dim(df_res)[2]
  L <- dim(df_res)[1]
  
  df_res <- as.matrix.data.frame(df_res)
  #df_res[,5:9] <- as.numeric(df_res[,5:9])
  #for(j in 1:5){
  # print(j)
  # df_res[,4+j] <- as.numeric(as.vector(df_res[,4+j]))
  #}  
  
  
  names_finder  <- function(dist, df_ene, dr) {
    
    tmp <- as.numeric(rownames(df_ene[df_ene$vet_dist < dist + dr & df_ene$vet_dist > dist - dr,]))
    
    choice <- sample(1:length(tmp),1, replace = T)
    
    chosen_couple <- tmp[choice]
    return(chosen_couple)
  }
  
  for(i in 1:N_links){
    print(i)
    
    #if( sum(df_res[,1]) != sum(df_res[,2])){break()}
    link <- sample(1:N_links, n, replace <- T)
    dist <- as.data.frame(df_ene$vet_dist[link])
    couples <- as.numeric(apply(dist, 1, names_finder, df_ene, dr))
    
    if(max(couples) > L || max(link) > L){break()}
    pos_1 <- c()
    pos_2 <- c()
    
    pos_1 <- (1:n)*L
    pos_2 <- pos_1
    pos_1 <- pos_1 + link
    pos_2 <- pos_2 + couples
    
    
    tmp_ene <- df_res[pos_1]
    df_res[pos_1] <- df_res[pos_2]
    df_res[pos_2] <- tmp_ene
    
    #df_res[link, 2:len] <- tmp_1
    #df_res[couples, 2:len] <- tmp_2
    
  }
  
  res <- as.data.frame(df_res[, 2:n+1])
  
  return(res)
  
}



##### 

strenth_local <- function(df_aus){
  
  EdgeList <- as.data.frame(df_aus)
  EdgeList[,1]=as.character(EdgeList[,1])
  EdgeList[,2]=as.character(EdgeList[,2])
  EdgeList <- as.matrix(EdgeList)
  g <- graph.edgelist(EdgeList[,1:2], directed=FALSE)
  E(g)$weight <- as.numeric(EdgeList[,3])
  
  Strength <- strength(g, weights = E(g)$weight)
  return(Strength)
}



###########################################################
###########################################################
##### Energy funtion 
###########################################################
###########################################################

energies <- function(pdb, LJcoeff, atomMap){
  
  Count_resNo <- 1
  
  for(i in 1:dim(pdb$atom)[1]){
    
    progress(i, max.value = dim(pdb$atom)[1], progress.bar = FALSE, char = "|",init = (i == 0), console = TRUE, gui = TRUE)
    
    
    ResNo <- as.numeric(gsub("[^0-9]", "", pdb$atom$resno[i]))
    
    if(i == 1){
      AnteResNo <- as.numeric(gsub("[^0-9]", "", pdb$atom$resno[i]))
      tmp <- gsub("[[:digit:]^]", "", pdb$atom$chain[i])
    }
    
    if(i != 1){
      AnteResNo <- as.numeric(gsub("[^0-9]", "", pdb$atom$resno[i-1]))
      tmp <- gsub("[[:digit:]^]", "", pdb$atom$chain[i-1])
    }
    
    if((ResNo < AnteResNo)|(as.character(pdb$atom$chain[i]) !=  tmp)){
      Count_resNo <- Count_resNo + 1 
    }
    
    pdb$atom$chain[i] <- paste0(pdb$atom$chain[i], Count_resNo)
    
    pdb$atom$elechain[i] <- paste(pdb$atom$resno[i], pdb$atom$insert[i], pdb$atom$resid[i], pdb$atom$chain[i], sep = "_")

        if(pdb$atom$elety[i] == "N" & pdb$atom$elety[i+1] == "HT1"){
      pdb$atom$resid[i] <- as.character("NTER")
    }
    if(pdb$atom$elety[i] == "C" & pdb$atom$elety[i+1] == "OT1"){
      pdb$atom$resid[i] <- as.character("CTER")
    }
    if(pdb$atom$elety[i] == "HT1" | pdb$atom$elety[i] == "HT2" | pdb$atom$elety[i] == "HT3"){
      pdb$atom$resid[i] <- as.character("NTER")
    }
    if(pdb$atom$elety[i] == "OT1" | pdb$atom$elety[i] == "OT2"){
      pdb$atom$resid[i] <- as.character("CTER")
    }
    for (j in 1:dim(atomMap)[1]) {
      if(pdb$atom$resid[i] == as.character(atomMap[j,"Resid"]) & pdb$atom$elety[i] == as.character(atomMap[j,"Atom"])){
        pdb$atom$newname[i] <- as.character(atomMap[j,"Name"])
        pdb$atom$charge[i] <- atomMap[j,"Charge"]
      }
    }
    
  }
  
  for (i in 1:dim(pdb$atom)[1]) {
    
    progress(i, max.value = dim(pdb$atom)[1], progress.bar = FALSE, char = "|",init = (i == 0), console = TRUE, gui = TRUE)
    
    
    for (j in 1:dim(LJcoeff)[1]) {
      if(pdb$atom$newname[i] == as.character(LJcoeff[j,"Names"])){
        pdb$atom$eps[i] <- LJcoeff[j,"eps"]
        pdb$atom$Rmin_2[i] <- as.character(LJcoeff[j,"Rmin_2"])
        pdb$atom$eps1_4[i] <- LJcoeff[j,"eps1_4"]
        pdb$atom$Rmin_2_1_4[i] <- LJcoeff[j,"Rmin_2_1_4"]
      }
      
    }
  } 
  
  Dist_Mat <- dm(pdb$xyz,grpby=pdb$atom[,'elechain'], scut = 2, ncore=1)
  rownames(Dist_Mat) <- unique(pdb$atom$elechain)  
  colnames(Dist_Mat) <- unique(pdb$atom$elechain)
  Dist_Mat <- Dist_Mat[rownames(Dist_Mat) %in% unique(pdb$atom$elechain[pdb$atom$type =="ATOM"]), colnames(Dist_Mat) %in% unique(pdb$atom$elechain[pdb$atom$type =="ATOM"])]
  
  dist_aa <- melt(Dist_Mat)
  dist_aa <- dist_aa[!is.na(dist_aa$value),]
  
  dist_aa <- dist_aa[dist_aa$value<=12,]
  
  colnames(dist_aa) <- c("X1", "X2", "value")
  
  vet_n <- c()
  vet_contacts <- c()
  vet_SumE <- c()
  vet_SumE_c <- c()
  vet_SumE_lj <- c()
  vet_dist <- c()
  
  for (n in unique(as.character(dist_aa$X1))){
    
    print(n)
    
    x_aus <- pdb$atom$x[pdb$atom$elechain==n]
    y_aus <- pdb$atom$y[pdb$atom$elechain==n]
    z_aus <- pdb$atom$z[pdb$atom$elechain==n]
    q_aus <- as.numeric(pdb$atom$charge[pdb$atom$elechain==n])
    eps_aus <- pdb$atom$eps[pdb$atom$elechain==n]
    Rmin_aus <- as.numeric(pdb$atom$Rmin_2[pdb$atom$elechain==n])
    
    contact <- as.character(dist_aa$X2[dist_aa$X1==n])
    
    if(length(contact)==0){
      
    }
    
    ene <- data.frame()
    
    for(i in 1:length(contact)){
      
      x <- pdb$atom$x[pdb$atom$elechain==contact[i]]
      y <- pdb$atom$y[pdb$atom$elechain==contact[i]]
      z <- pdb$atom$z[pdb$atom$elechain==contact[i]]
      q <- as.numeric(pdb$atom$charge[pdb$atom$elechain==contact[i]])
      eps <- pdb$atom$eps[pdb$atom$elechain==contact[i]]
      Rmin <- as.numeric(pdb$atom$Rmin_2[pdb$atom$elechain==contact[i]])
      E <- matrix(0,nrow = length(x_aus), ncol = length(x))
      E_c <- matrix(0,nrow = length(x_aus), ncol = length(x))
      E_lj <- matrix(0,nrow = length(x_aus), ncol = length(x))
      
      for(j in 1:length(x_aus)){
        
        for(k in 1:length(x)){
          
          r_c <- (sqrt((x_aus[j]-x[k])^2+(y_aus[j]-y[k])^2+(z_aus[j]-z[k])^2))*10^(-10)
          r <- sqrt((x_aus[j]-x[k])^2+(y_aus[j]-y[k])^2+(z_aus[j]-z[k])^2)
          cost <- 1/(4*pi*8.85418782*10^(-12))
          E_c[j,k] <-  (cost*(q_aus[j]*q[k]*2.56697*10^(-38))/r_c)*1.44*10^20
          E_lj[j,k] <- sqrt(eps_aus[j]*eps[k])*(((Rmin_aus[j]+Rmin[k])/r)^12-2*((Rmin_aus[j]+Rmin[k])/r)^6)
          E[j,k] <- E_c[j,k] + E_lj[j,k]
          
          
  
        }
      }
      
     
      
     
      
      vet_n <- c(vet_n, n)
      vet_contacts <- c(vet_contacts, contact[i])
      vet_SumE <- c(vet_SumE, sum(E))
      vet_SumE_c <-c(vet_SumE_c, sum(E_c))
      vet_SumE_lj <-c(vet_SumE_lj, sum(E_lj))
      vet_dist <- c(vet_dist, dist_aa[(dist_aa$X1==n | dist_aa$X2==n) &
                                        (dist_aa$X1 %in% contact[i] | dist_aa$X2 %in% contact[i]),3])
      ene <- data.frame(vet_n, vet_contacts, vet_dist, vet_SumE_c, vet_SumE_lj, vet_SumE)
      
    }
    
  }
  
  return(ene)
}



###########################################################
###########################################################
##### Energy funtionn atom - atom 
###########################################################
###########################################################


energies_atom_atom <- function(pdb_aus, LJcoeff, atomMap){
  
  pdb <- pdb_aus
  
  Count_resNo <- 1
  
  for (i in 1:dim(pdb$atom)[1]) {
    
    print(i)
    
    ResNo <- as.numeric(gsub("[^0-9]", "", pdb$atom$resno[i]))
    
    if(i == 1){
      AnteResNo <- as.numeric(gsub("[^0-9]", "", pdb$atom$resno[i]))
      tmp <- gsub("[[:digit:]^]", "", pdb$atom$chain[i])
    }
    
    if(i != 1){
      AnteResNo <- as.numeric(gsub("[^0-9]", "", pdb$atom$resno[i-1]))
      tmp <- gsub("[[:digit:]^]", "", pdb$atom$chain[i-1])
    }
    
    if((ResNo < AnteResNo)|(as.character(pdb$atom$chain[i]) !=  tmp)){
      Count_resNo <- Count_resNo + 1 
    }
    
    pdb$atom$chain[i] <- paste0(pdb$atom$chain[i], Count_resNo)
    
    pdb$atom$elechain[i] <- paste(pdb$atom$resno[i], pdb$atom$insert[i], pdb$atom$resid[i], pdb$atom$chain[i], pdb$atom$elety[i],sep = "_")
    
    if(pdb$atom$elety[i] == "N" & pdb$atom$elety[i+1] == "HT1"){
      pdb$atom$resid[i] <- as.character("NTER")
    }
    if(pdb$atom$elety[i] == "C" & pdb$atom$elety[i+1] == "OT1"){
      pdb$atom$resid[i] <- as.character("CTER")
    }
    if(pdb$atom$elety[i] == "HT1" | pdb$atom$elety[i] == "HT2" | pdb$atom$elety[i] == "HT3"){
      pdb$atom$resid[i] <- as.character("NTER")
    }
    if(pdb$atom$elety[i] == "OT1" | pdb$atom$elety[i] == "OT2"){
      pdb$atom$resid[i] <- as.character("CTER")
    }
    for (j in 1:dim(atomMap)[1]) {
      if(pdb$atom$resid[i] == as.character(atomMap[j,"Resid"]) & pdb$atom$elety[i] == as.character(atomMap[j,"Atom"])){
        pdb$atom$newname[i] <- as.character(atomMap[j,"Name"])
        pdb$atom$charge[i] <- atomMap[j,"Charge"]
      }
    }
    
  }
  
  for (i in 1:dim(pdb$atom)[1]) {
    
    print(i)
    
    for (j in 1:dim(LJcoeff)[1]) {
      if(pdb$atom$newname[i] == as.character(LJcoeff[j,"Names"])){
        pdb$atom$eps[i] <- LJcoeff[j,"eps"]
        pdb$atom$Rmin_2[i] <- as.character(LJcoeff[j,"Rmin_2"])
        pdb$atom$eps1_4[i] <- LJcoeff[j,"eps1_4"]
        pdb$atom$Rmin_2_1_4[i] <- LJcoeff[j,"Rmin_2_1_4"]
      }
      
    }
  } 
  
  Dist_Mat <- dm(pdb$xyz,grpby=pdb$atom[,'elechain'], scut = 2, ncore=6)
  rownames(Dist_Mat) <- unique(pdb$atom$elechain)  
  colnames(Dist_Mat) <- unique(pdb$atom$elechain)
  Dist_Mat <- Dist_Mat[rownames(Dist_Mat) %in% unique(pdb$atom$elechain[pdb$atom$type =="ATOM"]), colnames(Dist_Mat) %in% unique(pdb$atom$elechain[pdb$atom$type =="ATOM"])]
  
  dist_aa <- melt(Dist_Mat)
  dist_aa <- dist_aa[!is.na(dist_aa$value),]
  
  dist_aa <- dist_aa[dist_aa$value<=12,]
  
  ene_atom_atom <- function(row, pdb){ 
    #print("a")
    q_1 <- as.numeric(pdb$atom$charge[pdb$atom$elechain == row[1]])
    q_2<- as.numeric(pdb$atom$charge[pdb$atom$elechain == row[2]])
    
    
    eps_1 <- pdb$atom$eps[pdb$atom$elechain == row[1]]
    eps_2 <- pdb$atom$eps[pdb$atom$elechain == row[2]]
    
    Rm_1 <- as.numeric(pdb$atom$Rmin_2[pdb$atom$elechain == row[1]])  
    Rm_2 <- as.numeric(pdb$atom$Rmin_2[pdb$atom$elechain == row[2]])  
      
    r_c <- as.numeric(row[3])*10^(-10)
    r <- as.numeric(row[3])
    cost <- 1/(4*pi*8.85418782*10^(-12))
    E_c <-  (cost*(q_1*q_2*2.56697*10^(-38))/r_c)*1.44*10^20
    E_lj <- sqrt(eps_1*eps_2)*(((Rm_1+Rm_2)/r)^12-2*((Rm_1+Rm_2)/r)^6)
    E <- E_c + E_lj
  
    
    return(c(E_c,E_lj,E))
  
  }
  
  Energies <- apply(dist_aa, 1, ene_atom_atom,pdb )
  energia_atom <- as.data.frame(t(Energies))
  
  df_Energy_atom_atom <- as.data.frame(cbind(dist_aa, energia_atom))
  df_Energy_atom_atom_P1 <- (do.call("rbind", strsplit(as.character(df_Energy_atom_atom$Var1), "_")))
  df_Energy_atom_atom_P2 <- (do.call("rbind", strsplit(as.character(df_Energy_atom_atom$Var2), "_")))
  df_Energy_atom_atom_P3 <- (cbind(dist_aa$value, energia_atom))
  df_Energy_atom_atom_finale <- as.data.frame(cbind(df_Energy_atom_atom_P1,df_Energy_atom_atom_P2,df_Energy_atom_atom_P3))
  colnames(df_Energy_atom_atom_finale) <- c("Pos1","Insert1", "Res1","Chain1","Atom1","Pos2","Insert2", "Res2","Chain2","Atom2","Dist","E_C","E_LJ","E_Tot")
  
    
    ### conditions to remove rows linked to atom-atom interaction for the same or neighbor residues
    cond1 <- df_Energy_atom_atom_finale$Pos1 == df_Energy_atom_atom_finale$Pos2
    cond2 <- df_Energy_atom_atom_finale$Insert1 == df_Energy_atom_atom_finale$Insert2
    cond3 <- df_Energy_atom_atom_finale$Chain1 == df_Energy_atom_atom_finale$Chain2
  
    # 1?? condition
    df_return1 <- df_Energy_atom_atom_finale[ !(cond1 & cond2 & cond3),]
  
    cond1 <- df_return1$Pos1 == df_return1$Pos2
    cond2 <- df_return1$Insert1 == df_return1$Insert2
    cond3 <- df_return1$Chain1 == df_return1$Chain2
    cond4 <- as.numeric(as.vector(df_return1$Pos2)) - as.numeric(as.vector(df_return1$Pos1)) == 1
    cond5 <- as.numeric(as.vector(df_return1$Pos2)) - as.numeric(as.vector(df_return1$Pos1)) == -1
    
    # 2?? condition
    df_return2 <- df_return1[ ! ((cond4 & cond2 & cond3) | (cond5 & cond2 & cond3)) , ]
    
    
    cond6 <- df_return2$Insert1 == "" & df_return2$Insert2 == ""
    cond7 <- df_return2$Insert1 == "" & df_return2$Insert2 != ""
    cond8 <- df_return2$Insert1 != "" & df_return2$Insert2 == ""
    
    
    # 3?? condition
    df_clear <- df_return2[(cond6 | cond7 | cond8), ]
  
    
    # 4?? condition
    df_dirty <- df_return2[! (cond6 | cond7 | cond8), ]
    
    if (dim(df_dirty)[1] > 0){
      
   
    convert_letters_to_number <- function(vet_insert){
      
      diz <- cbind(LETTERS, 1:length(LETTERS))
      
      vet_num <- c()
      
      for(i in 1:length(vet_insert)){
        
        print(i)
        
        ele_aus <- vet_insert[i]
        num_aus <- diz[diz[,1] == ele_aus, 2]
        
        vet_num <- c(vet_num, num_aus)
      }
      
      return(vet_num)
    }
    
    
    df_dirty$Insert1_old <- df_dirty$Insert1
    df_dirty$Insert2_old <- df_dirty$Insert2
  
    df_dirty$Insert1 <- convert_letters_to_number(as.character(df_dirty$Insert1))
    df_dirty$Insert2 <- convert_letters_to_number(as.character(df_dirty$Insert2))
   
    cond9 <- abs(as.numeric(as.vector(df_dirty$Insert1)) - as.numeric(as.vector(df_dirty$Insert2))) == 1
    df_clear2 <- df_dirty[! cond9, ]
    
    df_return_final <- rbind(df_clear, df_clear2)
    
    return(df_return_final)
    }
    

  return(df_clear)
}









###########################################################
###########################################################
##### Network funtion 
###########################################################
###########################################################


Network_Generation <- function(Energy, Lower_Lim, Upper_Lim, Type){
  
  EdgeList <- as.data.frame(Energy[,c("vet_n", "vet_contacts", Type)])
  
  EdgeList <- EdgeList[(EdgeList[,colnames(EdgeList) == Type] < Upper_Lim) & (EdgeList[,colnames(EdgeList) == Type] > Lower_Lim), ]
  
  EdgeList[,1]=as.character(EdgeList[,1])
  EdgeList[,2]=as.character(EdgeList[,2])
  EdgeList <- as.matrix(EdgeList)
  g <- graph.edgelist(EdgeList[,1:2], directed=FALSE)
  E(g)$weight <- as.numeric(EdgeList[,3])
  
  return(g)
  
}


###########################################################
###########################################################
####### Local Network Parameters 
###########################################################
###########################################################


Local_Network_Paramters <- function(g, name_pdb){
 
  Betweenness_Centrality <- betweenness(g, v = V(g), directed = FALSE, nobigint = TRUE, normalized = TRUE,  weights = 1/E(g)$weight)
  Closeness_Centrality <- closeness(g, vids = V(g), weights = 1/E(g)$weight, normalized = TRUE)
  Strength_Value <- strength(g)
  Diversity_Index <- diversity(g, weights = E(g)$weight, vids = V(g))
  Mean_Shortest_Path <- apply(distances(g, algorithm = "dijkstra", weights = 1/E(g)$weight), 1, mean)
  Hub_Score_aus <- hub_score(g, weights = E(g)$weight)
  Hub_Score <- Hub_Score_aus$vector
  Clustering_Coefficient <- transitivity(g, type="barrat", weights = E(g)$weight)
  Degree_Value <- degree(g, loops = FALSE, normalized = TRUE)
  
  PDB_Name <- rep(name_pdb, length(Degree_Value))
  
  Parameters_DF <- cbind(PDB_Name, Betweenness_Centrality, Closeness_Centrality, Strength_Value, Diversity_Index, Mean_Shortest_Path, Hub_Score, Clustering_Coefficient, Degree_Value)
  Parameters_DF <- as.data.frame(Parameters_DF)
  
  return(Parameters_DF)
  
}



###########################################################
###########################################################
####### covert A in ALA and so on (and vice versa)
###########################################################
###########################################################

fn_conv_aa <- function(aa, condition){
  
  aa_3 <- c("ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HSD","ILE","LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL")
  aa_1 <- c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  
  df_aa_letters <- as.data.frame(cbind(aa_3,aa_1))
  colnames(df_aa_letters) <- c("aa_3","aa_1")
  
    if (condition == 0){
      return(as.character(df_aa_letters[df_aa_letters$aa_3==aa,2]))
    }
  
  if (condition == 1){
    return(as.character(df_aa_letters[df_aa_letters$aa_1==aa,1]))
  }
}




###########################################################
###########################################################
####### Assign the right type for each protein (Thermo, Meso)
###########################################################
###########################################################


Types_classificator <- function(dataset, Tresh_Meso, Tresh_Therm){
  
  dataset[dataset[,2] >= Tresh_Therm, 3] <- "T"
  dataset[dataset[,2] < Tresh_Meso, 3] <- "M"
  dataset[dataset[,2] > Tresh_Meso & dataset[,2] < Tresh_Therm, 3] <- "NC"
  colnames(dataset) <- c("PDB","Temp","Type")
  
  return(dataset)
  
}





###########################################################
###########################################################
####### Network Parameters Analysis
###########################################################
###########################################################

parameters_network_fn <- function(df_aus, lim_sx_1_neg, lim_dx_1_neg, lim_sx_1_pos, lim_dx_1_pos, lim_sx_2, lim_dx_2){
  
  ### Total Energy Network
  EdgeList <- as.data.frame(df_aus)
  EdgeList[,1]=as.character(EdgeList[,1])
  EdgeList[,2]=as.character(EdgeList[,2])
  EdgeList <- as.matrix(EdgeList)
  g <- graph.edgelist(EdgeList[,1:2], directed=FALSE)
  E(g)$weight <- as.numeric(EdgeList[,3])
  
  transitivity_g <- transitivity(g, type="barrat", weights = 1/E(g)$weight)
  degree_g <- degree(g, loops = FALSE, normalized = TRUE)
  vet_hub_score_g <- hub_score(g, weights = E(g)$weight)
  vet_hub_score_g <- vet_hub_score_g$vector
  vet_strength_g <- strength(g)
  
  df_Net_Par <- cbind(transitivity_g,degree_g,vet_hub_score_g,vet_strength_g, "Total")
  
  ### Strong Energy Network
  
  df_aus_strong <- df_aus[(df_aus$vet_SumE < lim_dx_1_pos & df_aus$vet_SumE > lim_sx_1_pos) | (df_aus$vet_SumE < lim_dx_1_neg & df_aus$vet_SumE > lim_sx_1_neg) | (df_aus$vet_SumE < lim_dx_2 & df_aus$vet_SumE > lim_sx_2),]
  
  EdgeList_strong <- as.data.frame(df_aus_strong)
  EdgeList_strong[,1]=as.character(EdgeList_strong[,1])
  EdgeList_strong[,2]=as.character(EdgeList_strong[,2])
  EdgeList_strong <- as.matrix(EdgeList_strong)
  g_strong <- graph.edgelist(EdgeList_strong[,1:2], directed=FALSE)
  E(g_strong)$weight <- as.numeric(EdgeList_strong[,3])
  
  transitivity_g_strong <- transitivity(g_strong, type="barrat", weights = 1/E(g_strong)$weight)
  degree_g_strong <- degree(g_strong, loops = FALSE, normalized = TRUE)
  vet_hub_score_g_strong <- hub_score(g_strong, weights = E(g_strong)$weight)
  vet_hub_score_g_strong <- vet_hub_score_g_strong$vector
  vet_strength_g_strong <- strength(g_strong)
  
  df_Net_Par_strong <- cbind(transitivity_g_strong,degree_g_strong,vet_hub_score_g_strong,vet_strength_g_strong, "Strong")
  
  
  ### Weak Energy Network
  
  df_aus_weak <- df_aus[(df_aus$vet_SumE < lim_sx_1_pos & df_aus$vet_SumE > lim_dx_1_neg),]
  
  EdgeList_weak <- as.data.frame(df_aus_weak)
  EdgeList_weak[,1]=as.character(EdgeList_weak[,1])
  EdgeList_weak[,2]=as.character(EdgeList_weak[,2])
  EdgeList_weak <- as.matrix(EdgeList_weak)
  g_weak <- graph.edgelist(EdgeList_weak[,1:2], directed=FALSE)
  E(g_weak)$weight <- as.numeric(EdgeList_weak[,3])
  
  transitivity_g_weak <- transitivity(g_weak, type="barrat", weights = 1/E(g_weak)$weight)
  degree_g_weak <- degree(g_weak, loops = FALSE, normalized = TRUE)
  vet_hub_score_g_weak <- hub_score(g_weak, weights = E(g_weak)$weight)
  vet_hub_score_g_weak <- vet_hub_score_g_weak$vector
  vet_strength_g_weak <- strength(g_weak)
  
  df_Net_Par_weak <- cbind(transitivity_g_weak,degree_g_weak,vet_hub_score_g_weak,vet_strength_g_weak, "Weak")
  
  
  ### output function
  df_Energy_Network <- rbind(df_Net_Par, df_Net_Par_strong, df_Net_Par_weak)
  
  return(df_Energy_Network)
}

##############
##### This function graphicates the D(S), random strength distribution showing left sdev and real S values
graficator <- function(df_rand){
  d <- dim(df_rand)[2]
  len <- dim(df_rand)[1]
  
  for (j in 2:d) {
    
    real <- mean(as.numeric(as.vector(df_rand[2:len,j])))
    sd_left <- real - sd(as.numeric(as.vector(df_rand[2:len,j])))
    
    min <- min(as.numeric(as.vector(df_rand[1:len, j])))
    max <- max(as.numeric(as.vector(df_rand[1:len, j])))
    
    min <- min - abs(1/10*min)
    max <- max + abs(1/10*max)
    
    hist(as.numeric(as.vector(df_rand[2:len,j])), 50,  main = paste(substr(FileList_rand_csv[i],1,4), colnames(df_rand)[j], sep = "---"), xlim=c(min,max), xlab="Strength", cex.lab=1.4, cex.axis=1.5)
    abline(v=df_rand[1,j], col="red", lwd=2)
    abline(v=real, col="black")
    abline(v=sd_left, col="black", lwd=2)
    
  }
}



##############
##### This function extracts the real Strength value, the mean and the standard deviation ones, from   the D(S), random strength distribution.
##############
gaussian_cumulant_computator <- function(df_rand){
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
  
  
  tmp <- as.data.frame(t(cbind(c(as.numeric(as.vector(vet_ave_true)),vet_ave_S, vet_sd_S))))
  
  real <- paste(colnames(df_rand[,2:d]), "real", sep = "_")
  mean <- paste(colnames(df_rand[,2:d]), "mean", sep = "_")
  sdev <- paste(colnames(df_rand[,2:d]), "sdev", sep = "_")
  
  
  colnames(tmp) <- c(real, mean, sdev)
  
  return(tmp)
}



###############
#### This function computes the ROC curve of a predictive method, given the predictor and the classification
###############

ROC_computator <- function(descriptor_vet, class_vet){
  
  
  pred <- prediction(descriptor_vet, class_vet)
  perf <- performance( pred , "sens", "spec")
  cut_pos<-which((perf@x.values[[1]]-perf@y.values[[1]])==min((perf@x.values[[1]]-perf@y.values[[1]])[perf@x.values[[1]]-perf@y.values[[1]]>0]))
  cutoff=unlist(pred@cutoffs)[cut_pos]
  plot(perf,colorize=TRUE,lty=2, lwd=2, cex.lab=1.4, cex.axis=1, las=2)
  auct <- paste(c("Cutoff = 1.0e-"),round(cutoff,1),sep="")
  legend(0.2,0.3,c(auct,"\n"),border="red",cex=0.7,box.col = "white")
  ####TP VS FP
  perf <- performance( pred , "tpr", "fpr")
  par(mar=c(5, 4, 4, 8) + 0.1)
  plot(perf,colorize=TRUE,lty=2, lwd=2, cex.lab=1.4, cex.axis=1, las=2)
  #mtext("-log10", side=4, line=3, cex.lab=0.5,las=2, col="blue")
  # calculating AUC
  auc <- round(performance(pred, "auc")@y.values[[1]],3)
  # now converting S4 class to vector
  auct <- paste(c("AUC = "),auc,sep="")
  legend(0.2,0.3,c(auct,"\n"),border="red",cex=1.5,box.col = "white")
  print(auc)
  return(auc)
}




####### This function computes the performance of the binary method given a threshold temperature, and the Roc curves

Stability_Analizator <- function(df_sample, Thres_Temp, ROC = 1){
  
  
  
  tmp <- Types_classificator(df_sample[, c("PDB", "Temp")], Thres_Temp,Thres_Temp)
  df_sample$Type <- tmp$Type
  
  
  ####### logical classification of proteins (useful for roc analysis etc...)
  
  df_sample[df_sample$Type=="T","Class"] <- 1
  df_sample[df_sample$Type=="M","Class"] <- 0
  
  
  dist_perc_1sd <- (df_sample[,1] - (df_sample[,2] - df_sample[,3]))/(df_sample[,2] - df_sample[,3])
  
  Roc <- NULL
  ###### ROC ANALYSIS #####
  if(ROC){   
    Roc <- ROC_computator(dist_perc_1sd, df_sample$Class)
  }
  
  
  
  for(i in 1:length(dist_perc_1sd)){
    if(dist_perc_1sd[i] > 0 & df_sample[i, "Type"] == 'T'){ df_sample[i, "Result"] <- 1}
    
    if(dist_perc_1sd[i] < 0 & df_sample[i, "Type"] == 'T'){ df_sample[i, "Result"] <- 0}
    
    if(dist_perc_1sd[i] > 0 & df_sample[i, "Type"] == 'M'){ df_sample[i, "Result"] <- 0}
    
    if(dist_perc_1sd[i] < 0 & df_sample[i, "Type"] == 'M'){ df_sample[i, "Result"] <- 1}
  }
  
  #### performance method
  performace <- dim(df_sample[df_sample[, "Result"] == 1,])[1]/dim(df_sample)[1]
  
  
  #plotting  Index vs Temp
  plot(dist_perc_1sd, df_sample$Temp, main = paste("Performace -> ", round(performace, 2),"\nRoc -> ", Roc,  sep = ''), xlim = c(min(dist_perc_1sd), max(dist_perc_1sd)), ylim=c(20,160), pch=20, ylab="Temperature", xlab="Theoretical Descriptor", cex=2, cex.axis=1.2, cex.lab=1.4)
  abline(v= 0, lwd=2, col = "red")
  abline(h = Thres_Temp, lwd=2, col="red")
  
  return(df_sample$Result)
}

