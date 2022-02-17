# BiPredict
Biclique based compound-protein-interaction (CPI) prediction


Bicliques consist of two types of nodes and edges connecting each node of different types (a). Here, blue circles represent compounds, c, and red squares represent proteins, p, while edges represent interactions. The biclique expansion starts with an existing biclique (b, yellow inner circle), here consisting of three compounds (c = 3) and two proteins (p = 2). Next, all compounds and proteins that are directly connected to any member of the biclique are identified. They represent interaction candidates (b, lightblue outer circle). Other compounds and proteins of the network which are not directly connected to any member of an existing biclique are not considered (b, grey squares and circles). Finally, interactions are predicted if interaction candidates lack only one edge to be a member of an existing biclique (b, green dashed lines). 
![methods_skizze](https://user-images.githubusercontent.com/82212543/126770821-6b673a8d-7bdc-4036-b748-0adf11509ded.png)

## Requirements
You need to install the following dependencies (R packages):
```
biclique
igraph
parallel
data.table 
```
The function 'bipredict' will compute all bicliques of the input CPI network using the R package 'biclique', predict
new interactions based on the biclique extension method and output a table of these new interactions.
If you want to get only predictions for bicliques above a minimum c/p-biclique size,
you can change the values of chemical_border and protein_border (otherwise c=4, p=2 are used).
The input file should be an edge list with two columns, one for compounds and one for proteins,
with columnnames "chemical" and "ProteinID" in the header of the input file (see example 'ecoli_edgelist.txt').
For large edge lists (>40,000 interactions), it is recommended to use the option 'Exclude_outliers = TRUE', which removes interactions to compounds or proteins which have a node degree of one, prior to biclique calculation. This step will reduce computational time, but at the cost that you will get a slightly lower number of predicted interactions for bicliques with only two nodes on any side of the biclique (c=2 or p=2). 

## Parameters
```
EDGE_LIST: Input edge list with two columns named: "chemical" and "ProteinID"
chemical_border: Least number of compounds (parameter c, explained above)
protein_border: Least number of proteins (parameter p, explained above)
CORES: Number of cores available for parallel computing (please check how many cores you have available, or set to 1, do not use all available cores)
Exclude_outliers: remove interactions to nodes with a node degree of one
```

## Usage 

We provide an additional script to compute an input edge list based on the KEGG database, which can directly be used as input for BiPredict. Please find it here: https://github.com/SandraThieme/KEGG_edge_list

```
source('bipredict.R')
edge_list = read.table('ecoli_edgelist.txt',header = T)
predicted_interactions_52 = bipredict(edge_list,chemical_border = 5,protein_border = 2,CORES = 1,Exclude_outliers=F)
write.table(predicted_interactions_52, file = 'my_predicted_interactions_52.txt', row.names = F)
```

Please find descriptions and examples also in the file bipredict.R You need to have both .R files in the same directory to run the predictions.

## Reference

Sandra Thieme and Dirk Walther.
Biclique extension as an effective approach to identify missing links in metabolic compound–protein interaction networks.
Bioinformatics Advances, Volume 2, Issue 1, 2022, vbac001
https://doi.org/10.1093/bioadv/vbac001

## Authors
Sandra Thieme
PhD Student
AG Bioinformatics

Max Planck Institute for Molecular Plant Physiology
Wissenschaftspark Golm, Am Mühlenberg 1, 14476 Potsdam-Golm 
https://www.mpimp-golm.mpg.de/7955/2walther
