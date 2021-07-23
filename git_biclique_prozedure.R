## ----------------------------------------------------------------
##
##   R package
##
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
##   02110-1301 USA
##
## -----------------------------------------------------------------

#' @title Biclique based CPI prediction
#'
#' @description
#' This function will compute the bicliques and output the statistics of these bicliques.
#' If you want to get bicliques above a threshold, you can change the values of lleast and rleast.
#' The input file should be tab delimited with number of vertices and edges at the head of the input file.
#' This package supports edge list.
#'
#' @param filename Input file name
#' @param left_least Least number of left partite <default = 1>
#' @param right_least Least number of right partite <default = 1>
#' @param version Algorithm version <default = 1>
#' @param filetype Input file format <default = 0>. 0-edge list, 1-binary matrix.
#' @param getclique Get bicliques <default = 1>. If you set it to 0. you'll only get the statistics without bicliques.
#' @param envir biclique environment
#'
#' @examples
#' bicliques = bi.clique(system.file("extdata", "example1.el", package = "biclique"))
#'
#' @export


library(biclique)
#https://cran.r-project.org/src/contrib/biclique_1.0.5.tar.gz
library(igraph)
#https://cran.r-project.org/src/contrib/igraph_1.2.6.tar.gz
library(data.table)
library(parallel)


source('git_calculate_bicliques.R',sep = '')


################################################################################
################################### options ####################################
################################################################################

# for parallel computing
CORES = 1
# save extended edgelist
#SAVE = F
# BiRewire random procedure
RANDOM = T
#real data without positives -> choose when reference positives list should be used
REAL = T
#EDGES = 1
Exclude_outliers = T
#select n percent positives from the network
P_POS = 5

###################################### biclique calculation ######################################
# the original interaction edge list from STITCH after selection unsing the confidence threshold
stitch_selected_edge_list = degree_edge_list_nonred_nondupl
rownames(stitch_selected_edge_list) <- NULL

if (Exclude_outliers){
  tab1 = stitch_selected_edge_list
  tab1$n = 1
  tab2 = aggregate(n ~ chemical,data = tab1,sum)
  tab3 = aggregate(n ~ ProteinID,data = tab1,sum)

  unique_chemicals = tab2[which(tab2$n == 1),"chemical"]
  unique_proteins = tab3[which(tab3$n == 1),"ProteinID"]
  unique_interactions = stitch_selected_edge_list[which(stitch_selected_edge_list$chemical %in% unique_chemicals | stitch_selected_edge_list$ProteinID %in% unique_proteins),]
  stitch_selected_edge_list = stitch_selected_edge_list[-as.numeric(rownames(unique_interactions)),]

  print(c('number of excluded unique proteins: ',length(unique_proteins)))
  print(c('number of excluded unique compounds: ',length(unique_chemicals)))
  print(c('number of interactions excluded: ',nrow(unique_interactions)))
  print(c('size of edgelist: ',nrow(stitch_selected_edge_list)))
  rm(tab1,tab2,tab3)
}
edge_list = stitch_selected_edge_list


#use the same validation edge list for each run
#select 10 percent edges which are used as positives and remove them for biclique calculation
percent_selected_edgelist = make_edgelist(stitch_selected_edge_list,percent = P_POS)
edge_list = percent_selected_edgelist[which(percent_selected_edgelist$selection == 0),c("chemical","ProteinID")]
rownames(edge_list)<-NULL
human_positives = percent_selected_edgelist[which(percent_selected_edgelist$selection == 1),c("chemical","ProteinID")]
rownames(human_positives)<-NULL

if(REAL){
  edge_list = stitch_selected_edge_list
}

# For RANDOM edgelist using BiRewire !!!
if (RANDOM){
  library(BiRewire)
  rnd_edge_list_1 = make_random_edgelist(edge_list)
  edge_list = rnd_edge_list_1
  colnames(edge_list) = c('chemical','ProteinID')
}

# you can specify (EDGE_LIST,left_border,right_border)
# default values are (left_border = 2,right_border = 2)
print('calculate bicliques')
my_bicliques = calculate_bicliques(CORES,edge_list)
my_bicliques = my_bicliques[!is.na(my_bicliques)]
print(timestamp())
all_my_bicliques = list()
for (i in 1:length(my_bicliques)){
  all_my_bicliques = append(all_my_bicliques,my_bicliques[[i]])
}
names(all_my_bicliques) = paste(seq(1,length(all_my_bicliques)),names(all_my_bicliques),sep = '')
if (is.na(all_my_bicliques[1])){
  print('No bicliques are contained in the input data.')
  stop()
}else{
  print('extend bicliques')
  print(timestamp())
  if (EDGES == 1){
    extended_bicliques = extend_biclique(all_my_bicliques,edge_list,CORES)
    rm(all_my_bicliques)
  }
  else if(EDGES == 2){
    extended_bicliques = extend_biclique_2(all_my_bicliques,edge_list,CORES)
    rm(all_my_bicliques)
  }
  else {print('Edges attribute not specified.')}
  print(timestamp())
}

tx = timestamp()
extended_bicliques = main(edge_list,CORES)
if (REAL){
  #make result tab for real data
  edge_list_extended = edgelist_from_bicliques(extended_bicliques,edge_list,CORES)
  print('make result tab for real data without positives dataset.')
  result_tab_real_human = make_result_tab_real(CORES,edge_list,extended_bicliques,human_negatives,max_size = 10)
}else{
  #results tab with positive data set removed from the edge list
  edge_list_extended = edgelist_from_bicliques(extended_bicliques,edge_list,CORES)
  result_tab = make_result_tab_parallel(edge_list_extended,CORES,edge_list,extended_bicliques,human_positives,human_negatives,max_size = 10)
}

#print end time
print(tx)
timestamp()
rm(tx)


