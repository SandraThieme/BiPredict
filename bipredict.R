source('extend_bicliques.R')

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
#' This function will compute all bicliques of the input CPI network using the R package 'biclique', predict
#' new interactions based on the biclique extension method and output a table of these new interactions.
#' If you want to get only predictions for bicliques above a munimum c/p-biclique size,
#' you can change the values of chemical_border and protein_border.
#' The input file should be a edge list with two columns, one for compounds and one for proteis,
#' with columnames "chemical" and "ProteinID" at the head of the input file. (see example 'ecoli_edgelist.txt')
#'
#' @param EDGE_LIST Input edge list
#' @param chemical_border Least number of compounds
#' @param protein_border Least number of proteins
#' @param CORES number of used cores for parallel computing (please check how much cores you have available, or use 1, never use all available cores)
#' @examples
#' edge_list = read.table('ecoli_edgelist.txt',header = T)
#' predicted_interactions_52 = bipredict(edge_list,5,2,1)


library(biclique)
#https://cran.r-project.org/src/contrib/biclique_1.0.5.tar.gz
library(igraph)
#https://cran.r-project.org/src/contrib/igraph_1.2.6.tar.gz
library(data.table)
library(parallel)


bipredict = function (EDGE_LIST,chemical_border = 4,protein_border = 2,CORES=1){
  my_bicliques = calculate_bicliques(CORES,EDGE_LIST,chemical_border,protein_border)
  my_extended_bicliques = extend_biclique(my_bicliques,EDGE_LIST,CORES)
  predicted_interactions = edgelist_from_bicliques(my_extended_bicliques,EDGE_LIST,CORES)
  return(predicted_interactions)
}
