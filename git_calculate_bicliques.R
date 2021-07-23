library(biclique)
library(data.table)
library(parallel)
library(igraph)

# make edge list with randomly selected x percent as positives
make_edgelist = function(stitch_selected_edge_list,percent = 10){
  rownames(stitch_selected_edge_list)<-NULL
  colnames(stitch_selected_edge_list) = c('chemical','ProteinID')
  vec10 = as.numeric(sample(rownames(stitch_selected_edge_list),round(nrow(stitch_selected_edge_list)/100*percent)))
  percent_selected_edgelist = stitch_selected_edge_list
  percent_selected_edgelist$selection = 0
  percent_selected_edgelist[vec10,'selection'] = 1
  return(percent_selected_edgelist)
}

# calculate bicliques using R package biclique, parallel if CORES > 1 for each connected component of the input network
calculate_bicliques = function (CORES,EDGE_LIST,chemical_border = 2,protein_border = 2){
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
  return(my_bicliques)
}

# extend bicliques = predict candidate interactions, parallel for the bicliques if CORES > 1
extend_biclique = function(my_bicliques,working_edge_list,CORES){
  colnames(working_edge_list) = c('chemical','ProteinID')
  print('start')

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

# use this function?
extend_biclique_list = function(my_bicliques,working_edge_list,CORES){
  colnames(working_edge_list) = c('chemical','ProteinID')
  print('start')

  # function for parallel interaction candidate calculation
  do_it = function(X){
    one_biclique = my_bicliques
    clique_one = c(one_biclique$left,one_biclique$right)
    one_biclique$clique_one = clique_one
    #print('step1')
    biclique_edgelist = working_edge_list[which(working_edge_list$chemical %in% clique_one | working_edge_list$ProteinID %in% clique_one),]
    biclique_edgelist = unique(biclique_edgelist[,c('chemical','ProteinID')])
    biclique_edgelist$n = 1
    biclique_edgelist_chemicals = aggregate(n ~ chemical,biclique_edgelist,sum)
    biclique_edgelist_proteins = aggregate(n ~ ProteinID,biclique_edgelist,sum)
    #print('step2')
    left_min_clique_length = length(one_biclique$left)-1
    right_min_clique_length = length(one_biclique$right)-1
    left_candidates = setdiff(biclique_edgelist_chemicals[which(biclique_edgelist_chemicals$n >=right_min_clique_length),"chemical"],clique_one)
    right_candidates = setdiff(biclique_edgelist_proteins[which(biclique_edgelist_proteins$n >=left_min_clique_length),"ProteinID"],clique_one)
    one_biclique$left_candidates = left_candidates
    one_biclique$right_candidates = right_candidates
    return(one_biclique)
  }
  my_bicliques = lapply(X = my_bicliques,FUN = do_it)
  return(my_bicliques)
}

#make edge list for extended bicliques
edgelist_from_bicliques = function(extended_biclique,edge_list,CORES){
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
  colnames(edge_list) = colnames(my_edge_list)[1:2]

  library(data.table)
  my_edge_list = data.table(my_edge_list,key = colnames(my_edge_list[1:2]))
  edge_list = data.table(edge_list,key = colnames(edge_list[1:2]))
  edge_list$Var3 = 'stitch'
  merged_edgelist = unique(merge(edge_list,my_edge_list,all.y = T))

  # edge_list$Var3 = 'stitch'
  # merged_edgelist = unique(merge(edge_list,my_edge_list,all.y = T))

  merged_edgelist = merged_edgelist[is.na(merged_edgelist$Var3),]
  # edge_list$Var3<-NULL
  merged_edgelist$Var3<-NULL
  my_edge_list = as.data.frame(merged_edgelist)
  print(nrow(my_edge_list))

  #return(unique(my_edge_list))
  return(my_edge_list)
}

# table of performance measures
make_result_tab = function(CORES,edge_list,extended_bicliques,edge_list_positive,known_edge_list_negative,max_size = 10){
  colnames(edge_list) = c('chemical','ProteinID')
  working_edgelist = edge_list

  ### faster if you have many cores for parallel computing ###
  edge_list_extended = edgelist_from_bicliques(extended_bicliques,edge_list,CORES)
  colnames(edge_list_extended) = c('chemical','ProteinID','biclique')

  result_tab = data.frame()

  n_compounds = lapply(extended_bicliques, FUN = function(x) length(x[['left']]))
  n_proteins = lapply(extended_bicliques, FUN = function(x) length(x[['right']]))
  biclique_size = data.frame(seq(1:length(extended_bicliques)),unlist(n_compounds),unlist(n_proteins))
  colnames(biclique_size) = c('number','compounds','proteins')
  biclique_size$compounds = as.numeric(biclique_size$compounds)
  biclique_size$proteins = as.numeric(biclique_size$proteins)
  biclique_size$number = as.numeric(biclique_size$number)
  biclique_size$size = biclique_size$compounds*biclique_size$proteins

  left_border  = c(9,9,8,9,8,8,8,8,7,7,7,6,6,6,6,6,5,5,5,5,4,4,4,4,4,3,3,3,2,2)
  right_border = c(6,5,9,8,4,5,3,2,4,3,2,2,3,4,5,6,2,3,4,5,3,2,4,5,7,9,2,3,2,9)

  # apply only biclique sizes which are contained in the actual data
  all_biclique_sizes = unique(biclique_size[,c('compounds','proteins')])
  all_biclique_sizes = all_biclique_sizes[which(all_biclique_sizes$compounds < max_size & all_biclique_sizes$proteins < max_size),]
  rownames(all_biclique_sizes)<-NULL

  # all sizes
  # left_border = as.vector(all_biclique_sizes$compounds)
  # right_border = as.vector(all_biclique_sizes$proteins)

  # only selection of sizes
  biclique_size_borders = cbind(left_border,right_border)
  colnames(biclique_size_borders) = colnames(all_biclique_sizes)
  biclique_size_borders = unique(merge(all_biclique_sizes,biclique_size_borders))
  left_border = as.vector(biclique_size_borders$compounds)
  right_border = as.vector(biclique_size_borders$proteins)

  j = 1
  print(length(left_border))
  for (i in 1:length(left_border)){
    print(i)
    # make edge list only including bicliques of the selected sizes or of bigger size (more edges), out of extended bicliques list
    edge_number = left_border[i]*right_border[i]
    #biclique_vec = biclique_size[which(biclique_size$size >= edge_number),"number"]
    biclique_vec = biclique_size[which(biclique_size$compounds >= left_border[i] & biclique_size$proteins >= right_border[i]),"number"]


    # edge list of predicted interactions
    #selected_edge_list_extended = edgelist_from_biclique_list(extended_bicliques[biclique_vec],edge_list)

    selected_edge_list_extended = edge_list_extended[which(edge_list_extended$biclique %in% biclique_vec),c('chemical','ProteinID')]
    selected_extended_bicliques = selected_edge_list_extended
    if(nrow(selected_edge_list_extended)==0){print(c('No predictions for biclique size: ',left_border[i],right_border[i]));next}

    selected_edge_list_extended = unique(selected_edge_list_extended)
    selected_edge_list_extended$n = 1
    selected_edge_list_extended = aggregate(n ~. , data = selected_edge_list_extended,sum)
    # print(c('Number of predicted edges: ',sum(selected_edge_list_extended$n)))
    # predictions = sum(selected_edge_list_extended$n) # sum of all predictions, including known and multiple predicted edges
    # pred_edges = nrow(selected_edge_list_extended) # unique number of all predictions

    #print(c('Number of predicted edges without previously known edges: ',sum(selected_edge_list_extended$n)))
    #predictions_unknown = sum(selected_edge_list_extended$n) # sum of predictions, without known edges including multiple predictions
    pred_edges_unknown = nrow(selected_edge_list_extended) # unique number of unknown predictions, previously known edges were removed
    selected_edge_list_extended$n<-NULL

    #positives which are contained in original edge list (relevant for paper data - not realy for random sampled positives)
    tab1 = edge_list_positive[,c("chemical","ProteinID")]
    colnames(tab1) = colnames(edge_list)
    tab_pos = unique(merge(edge_list,tab1))

    tab1 = known_edge_list_negative[,c("chemical","ProteinID")]
    colnames(tab1) = colnames(edge_list)
    tab_neg = unique(merge(edge_list,tab1))

    tab1 = edge_list_positive[,c("chemical","ProteinID")]
    colnames(tab1) = colnames(selected_edge_list_extended)
    tab_ext_pos = unique(merge(selected_edge_list_extended,tab1))

    tab1 = known_edge_list_negative[,c("chemical","ProteinID")]
    colnames(tab1) = colnames(selected_edge_list_extended)
    tab_ext_neg = unique(merge(selected_edge_list_extended,tab1))

    # number of positives which were possible to be predicted
    px = nrow(edge_list_positive[which(edge_list_positive$ProteinID %in% working_edgelist$ProteinID & edge_list_positive$chemical %in% working_edgelist$chemical),])
    # number of negatives which were possible to be predicted
    py = nrow(known_edge_list_negative[which(known_edge_list_negative$ProteinID %in% working_edgelist$ProteinID & known_edge_list_negative$chemical %in% working_edgelist$chemical),])

    result_tab[j,'Dataset'] = paste('Biclique c:',left_border[i],', p:',right_border[i],sep = '')
    result_tab[j,'Biclique.size.min'] = edge_number
    result_tab[j,'Bicliques'] = length(biclique_vec)
    # result_tab[j,'Predictions'] = predictions
    # result_tab[j,'Predicted.edges'] = pred_edges
    # result_tab[j,'Novel.predictions'] = predictions_unknown
    result_tab[j,'Edges'] = pred_edges_unknown
    result_tab[j,'px'] = px - nrow(tab_pos)
    result_tab[j,'TP'] = nrow(tab_ext_pos)
    result_tab[j,'TPR'] = round(nrow(tab_ext_pos)/((px-nrow(tab_pos))/100),2)
    result_tab[j,'py'] = py - nrow(tab_neg)
    result_tab[j,'FP'] = nrow(tab_ext_neg)
    result_tab[j,'FPR'] = round(nrow(tab_ext_neg)/((py-nrow(tab_neg))/100),2)
    result_tab[j,'PPV'] = round((result_tab[j,'TP']/(result_tab[j,'TP']+result_tab[j,'FP']))*100,4)
    result_tab = result_tab[order(result_tab$Edges),]
    rownames(result_tab)<-NULL

    # counter for rows in result table
    j = j + 1
  }
  return(result_tab)
}

# table of performance measures
make_result_tab_real = function(CORES,edge_list,extended_bicliques,known_edge_list_negative,max_size = 10){
  working_edgelist = edge_list
  colnames(working_edgelist) = c('chemical','ProteinID')
  edge_list_extended = edgelist_from_bicliques(extended_bicliques,edge_list,CORES)
  colnames(edge_list_extended) = c('chemical','ProteinID','biclique')

  result_tab = data.frame()
  biclique_size = data.frame()
  #left are chemicals, right are proteins
  biclique_size[1:length(extended_bicliques),'number'] = seq(1,length(extended_bicliques))
  for (i in 1:length(extended_bicliques)){
    biclique_size[i,'left'] = length(extended_bicliques[[i]]$left)
    biclique_size[i,'right'] = length(extended_bicliques[[i]]$right)
  }

  # applied minumum biclique sizes for chemicals and proteins
  # left_border = c(6,6,6,6,6,5,5,5,5,4,4,4,3,3,2)
  # right_border = c(2,3,4,5,6,2,3,4,5,2,3,4,2,3,2)
  #reduced when using no degree cutoff
  left_border  = c(8,8,8,7,7,7,6,6,6,6,6,5,5,5,5,4,4,4,3,3,2)
  right_border = c(4,3,2,4,3,2,2,3,4,5,6,2,3,4,5,3,2,4,2,3,2)

  # apply only bilcique sizes which are contained in the actual data
  all_biclique_sizes = unique(biclique_size[,c("left","right")])
  all_biclique_sizes = all_biclique_sizes[which(all_biclique_sizes$left <= max_size & all_biclique_sizes$right <= max_size),]
  biclique_size_borders = cbind(left_border,right_border)
  colnames(biclique_size_borders) = colnames(all_biclique_sizes)
  biclique_size_borders = unique(merge(all_biclique_sizes,biclique_size_borders))
  left_border = as.vector(biclique_size_borders$left)
  right_border = as.vector(biclique_size_borders$right)

  j = 1
  print(length(left_border))
  for (i in 1:length(left_border)){
    print(i)
    # make edge list only including bicliques of the selected minimum sizes, out of extended bicliques list
    biclique_vec = biclique_size[which(biclique_size$left >= left_border[i] & biclique_size$right >= right_border[i]),"number"]
    selected_edge_list_extended = edge_list_extended[which(edge_list_extended$biclique %in% biclique_vec),c('chemical','ProteinID')]
    selected_extended_bicliques = selected_edge_list_extended
    if(nrow(selected_extended_bicliques)==0){print(c('No predictions for biclique size: ',left_border[i],right_border[i]));next}

    #selected_edge_list_extended = unique(selected_edge_list_extended)
    selected_edge_list_extended$n = 1
    selected_edge_list_extended = aggregate(n ~. , data = selected_edge_list_extended,sum)
    #print(c('Number of predicted edges: ',sum(selected_edge_list_extended$n)))
    # predictions = sum(selected_edge_list_extended$n) # sum of all predictions, including known and multiple predicted edges
    # pred_edges = nrow(selected_edge_list_extended) # unique number of all predictions

    #print(c('Number of predicted edges without previously known edges: ',sum(selected_edge_list_extended$n)))
    predictions_unknown = sum(selected_edge_list_extended$n) # sum of predictions, without known edges including multiple predictions
    pred_edges_unknown = nrow(selected_edge_list_extended) # unique number of unknown predictions, previously known edges were removed
    selected_edge_list_extended$n<-NULL

    tab1 = known_edge_list_negative[,c("chemical","ProteinID")]
    colnames(tab1) = colnames(edge_list)
    tab_neg = unique(merge(edge_list,tab1))

    tab1 = known_edge_list_negative[,c("chemical","ProteinID")]
    colnames(tab1) = colnames(selected_edge_list_extended)
    tab_ext_neg = unique(merge(selected_edge_list_extended,tab1))

    # number of ngatives which were possible to be predicted
    py = nrow(known_edge_list_negative[which(known_edge_list_negative$ProteinID %in% working_edgelist$ProteinID & known_edge_list_negative$chemical %in% working_edgelist$chemical),])
    result_tab[j,'Dataset'] = paste('Biclique c:',left_border[i],', p:',right_border[i],sep = '')
    result_tab[j,'Bicliques'] = length(biclique_vec)
    # result_tab[j,'Predictions'] = predictions
    # result_tab[j,'Predicted.edges'] = pred_edges
    result_tab[j,'Novel.predictions'] = predictions_unknown
    result_tab[j,'Edges'] = pred_edges_unknown
    result_tab[j,'py'] = py-nrow(tab_neg)
    result_tab[j,'Negatives'] = nrow(tab_ext_neg)
    result_tab[j,'FPR'] = round(nrow(tab_ext_neg)/((py-nrow(tab_neg))/100),2)
    result_tab = result_tab[order(result_tab$FPR),]
    rownames(result_tab)<-NULL
    # counter for rows in result table
    j = j + 1
  }
  return(result_tab)
}

# table of performance measures
make_result_tab_parallel = function(edge_list_extended,CORES,edge_list,extended_bicliques,edge_list_positive,known_edge_list_negative,max_size = 10){
  library(parallel)
  colnames(edge_list) = c('chemical','ProteinID')
  working_edgelist = edge_list
  colnames(edge_list_extended) = c('chemical','ProteinID','biclique')
  result_tab = data.frame()

  n_compounds = lapply(extended_bicliques, FUN = function(x) length(x[['left']]))
  n_proteins = lapply(extended_bicliques, FUN = function(x) length(x[['right']]))
  biclique_size = data.frame(seq(1:length(extended_bicliques)),unlist(n_compounds),unlist(n_proteins))
  colnames(biclique_size) = c('number','compounds','proteins')
  biclique_size$compounds = as.numeric(biclique_size$compounds)
  biclique_size$proteins = as.numeric(biclique_size$proteins)
  biclique_size$number = as.numeric(biclique_size$number)
  biclique_size$size = biclique_size$compounds*biclique_size$proteins

  # applied minimum biclique sizes for chemicals and proteins
  left_border  = c(9,9,8,9,8,8,8,8,7,7,7,6,6,6,6,6,5,5,5,5,4,4,4,4,4,3,3,3,2,2)
  right_border = c(6,5,9,8,4,5,3,2,4,3,2,2,3,4,5,6,2,3,4,5,3,2,4,5,7,9,2,3,2,9)
  # all sizes
  # left_border = as.vector(all_biclique_sizes$compounds)
  # right_border = as.vector(all_biclique_sizes$proteins)

  # apply only biclique sizes which are contained in the data
  all_biclique_sizes = unique(biclique_size[,c('compounds','proteins')])
  all_biclique_sizes = all_biclique_sizes[which(all_biclique_sizes$compounds < max_size & all_biclique_sizes$proteins < max_size),]
  rownames(all_biclique_sizes)<-NULL

  # only selection of sizes
  biclique_size_borders = cbind(left_border,right_border)
  colnames(biclique_size_borders) = colnames(all_biclique_sizes)
  biclique_size_borders = unique(merge(all_biclique_sizes,biclique_size_borders))
  left_border = as.vector(biclique_size_borders$compounds)
  right_border = as.vector(biclique_size_borders$proteins)

  j = 1
  print(length(left_border))

  do_it = function(i){
    print(i)
    # make edge list only including bicliques of the selected sizes or of bigger size (more edges), out of extended bicliques list
    edge_number = left_border[i]*right_border[i]
    biclique_vec = biclique_size[which(biclique_size$compounds >= left_border[i] & biclique_size$proteins >= right_border[i]),"number"]

    selected_edge_list_extended = edge_list_extended[which(edge_list_extended$biclique %in% biclique_vec),c('chemical','ProteinID')]
    selected_extended_bicliques = selected_edge_list_extended
    if(nrow(selected_edge_list_extended)==0){print(c('No predictions for biclique size: ',left_border[i],right_border[i]));next}

    selected_edge_list_extended = unique(selected_edge_list_extended)
    selected_edge_list_extended$n = 1
    selected_edge_list_extended = aggregate(n ~. , data = selected_edge_list_extended,sum)
    pred_edges_unknown = nrow(selected_edge_list_extended) # unique number of unknown predictions, previously known edges were removed
    selected_edge_list_extended$n<-NULL

    #positives which are contained in original edge list
    tab1 = edge_list_positive[,c("chemical","ProteinID")]
    colnames(tab1) = colnames(edge_list)
    tab_pos = unique(merge(edge_list,tab1))

    tab1 = known_edge_list_negative[,c("chemical","ProteinID")]
    colnames(tab1) = colnames(edge_list)
    tab_neg = unique(merge(edge_list,tab1))

    tab1 = edge_list_positive[,c("chemical","ProteinID")]
    colnames(tab1) = colnames(selected_edge_list_extended)
    tab_ext_pos = unique(merge(selected_edge_list_extended,tab1))

    tab1 = known_edge_list_negative[,c("chemical","ProteinID")]
    colnames(tab1) = colnames(selected_edge_list_extended)
    tab_ext_neg = unique(merge(selected_edge_list_extended,tab1))

    # number of positives which were possible to be predicted
    px = nrow(edge_list_positive[which(edge_list_positive$ProteinID %in% working_edgelist$ProteinID & edge_list_positive$chemical %in% working_edgelist$chemical),])
    # number of negatives which were possible to be predicted
    py = nrow(known_edge_list_negative[which(known_edge_list_negative$ProteinID %in% working_edgelist$ProteinID & known_edge_list_negative$chemical %in% working_edgelist$chemical),])

    result_tab[j,'Dataset'] = paste('Biclique c:',left_border[i],', p:',right_border[i],sep = '')
    result_tab[j,'Biclique.size.min'] = edge_number
    result_tab[j,'Bicliques'] = length(biclique_vec)
    result_tab[j,'Edges'] = pred_edges_unknown
    result_tab[j,'px'] = px - nrow(tab_pos)
    result_tab[j,'TP'] = nrow(tab_ext_pos)
    result_tab[j,'TPR'] = round(nrow(tab_ext_pos)/((px-nrow(tab_pos))/100),2)
    result_tab[j,'py'] = py - nrow(tab_neg)
    result_tab[j,'FP'] = nrow(tab_ext_neg)
    result_tab[j,'FPR'] = round(nrow(tab_ext_neg)/((py-nrow(tab_neg))/100),2)
    result_tab[j,'PPV'] = round((result_tab[j,'TP']/(result_tab[j,'TP']+result_tab[j,'FP']))*100,4)
    result_tab = result_tab[order(result_tab$Edges),]
    rownames(result_tab)<-NULL

    return(result_tab)
  }
  result_tab_list = mclapply(seq(1:length(left_border)),FUN = do_it,mc.cores = CORES)
  result_tab = do.call(rbind,result_tab_list)

  return(result_tab)
}

