library(biclique)
library(data.table)
library(parallel)
library(igraph)


# calculate bicliques using R package biclique, parallel if CORES > 1 for each connected component of the input network
calculate_bicliques = function (CORES,EDGE_LIST,chemical_border = 2,protein_border = 2){
  print('calculate bicliques')
  edge_list = EDGE_LIST
  edge_matrix = unique(as.matrix(edge_list[,c("chemical","ProteinID")]))
  rownames(edge_matrix)<-NULL
  graph2 = graph_from_edgelist(edge_matrix,directed = F)
  graph2 = igraph::simplify(graph2)
  graph2_components = components(graph2)

  # consider only components with more than three members
  x = names(table(graph2_components$membership)[which(table(graph2_components$membership)>3)])
  print(length(x))
  # function for parallel biclique calculation for each connected component
  do_it = function(i){
    PATH_TO_EDGELIST = tempfile(pattern = "", fileext = paste('_',i,".el",sep = ''))
    x_ids = names(which(graph2_components$membership == i))
    EDGE_LIST = edge_list[which(edge_list$chemical %in% x_ids & edge_list$ProteinID %in% x_ids),]
    write.table(EDGE_LIST,quote = F,row.names = F,col.names = F,file = PATH_TO_EDGELIST,sep = '\t')
    bi.format(PATH_TO_EDGELIST)
    my_bicliques <- tryCatch(
      {bi.clique(PATH_TO_EDGELIST,left_least = chemical_border,right_least = protein_border)}, error=function(cond){return(NA)})
    return(my_bicliques)
  }
  my_bicliques = mclapply(X = x, FUN = do_it, mc.cores = CORES)

  my_bicliques = my_bicliques[!is.na(my_bicliques)]
  all_my_bicliques = list()
  for (i in 1:length(my_bicliques)){
    all_my_bicliques = append(all_my_bicliques,my_bicliques[[i]])
  }
  names(all_my_bicliques) = paste(seq(1,length(all_my_bicliques)),names(all_my_bicliques),sep = '')
  if (is.na(all_my_bicliques[1])){
    print('No bicliques are contained in the input data.')
    stop()
  }

  return(all_my_bicliques)
}

# extend bicliques = predict candidate interactions, parallel for the bicliques if CORES > 1
extend_biclique = function(my_bicliques,working_edge_list,CORES){
  colnames(working_edge_list) = c('chemical','ProteinID')
  print('extend bicliques')

  # function for parallel interaction candidate calculation
  do_it = function(name){
    one_biclique = my_bicliques[name]
    clique_one = c(one_biclique[[name]]$left,one_biclique[[name]]$right)
    one_biclique[[name]]$clique_one = clique_one
    #print('step1')
    biclique_edgelist = working_edge_list[which(working_edge_list$chemical %in% clique_one | working_edge_list$ProteinID %in% clique_one),]
    biclique_edgelist = unique(biclique_edgelist[,c('chemical','ProteinID')])
    biclique_edgelist$n = 1
    biclique_edgelist_chemicals = aggregate(n ~ chemical,biclique_edgelist,sum)
    biclique_edgelist_proteins = aggregate(n ~ ProteinID,biclique_edgelist,sum)
    #print('step2')
    left_min_clique_length = length(one_biclique[[name]]$left)-1
    right_min_clique_length = length(one_biclique[[name]]$right)-1
    left_candidates = setdiff(biclique_edgelist_chemicals[which(biclique_edgelist_chemicals$n >=right_min_clique_length),"chemical"],clique_one)
    right_candidates = setdiff(biclique_edgelist_proteins[which(biclique_edgelist_proteins$n >=left_min_clique_length),"ProteinID"],clique_one)
    one_biclique[[name]]$left_candidates = left_candidates
    one_biclique[[name]]$right_candidates = right_candidates
    return(one_biclique[[name]])
  }
  my_bicliques = mclapply(X = names(my_bicliques),FUN = do_it,mc.cores = CORES)
  return(my_bicliques)
}

#make edge list for extended bicliques
edgelist_from_bicliques = function(extended_biclique,edge_list,CORES){
  print('make edgelist')
  library(parallel)
  do_it = function(i){
    temp_edge_list_one = expand.grid(extended_biclique[[i]]$left, extended_biclique[[i]]$right_candidates,stringsAsFactors =F)
    temp_edge_list_two = expand.grid(extended_biclique[[i]]$left_candidates, extended_biclique[[i]]$right,stringsAsFactors =F)
    temp_edge_list = rbind(temp_edge_list_one,temp_edge_list_two)
    if(nrow(temp_edge_list)>0){temp_edge_list$biclique = i}
    return(temp_edge_list)
  }
  my_edge_list = do.call(rbind,mclapply(X = seq(1:length(extended_biclique)),FUN = do_it,mc.cores = CORES))

  # delete entries in extended edgelist which are already in original edgelist (previously known edges)
  print(nrow(my_edge_list))
  my_edge_list = unique(my_edge_list[,c('Var1','Var2')])
  colnames(my_edge_list) = colnames(edge_list)[1:2]

  library(data.table)
  my_edge_list = data.table(my_edge_list,key = colnames(my_edge_list[1:2]))
  edge_list = data.table(edge_list,key = colnames(edge_list[1:2]))
  edge_list$Var3 = 'stitch'
  merged_edgelist = unique(merge(edge_list,my_edge_list,all.y = T))

  merged_edgelist = merged_edgelist[is.na(merged_edgelist$Var3),]
  merged_edgelist$Var3<-NULL
  my_edge_list = as.data.frame(merged_edgelist)
  print(nrow(my_edge_list))
  #colnames(my_edge_list) = c('chemical','ProteinID','biclique')
  #my_edge_list = unique(my_edge_list[,c('chemical','ProteinID')])

  return(my_edge_list)
}
