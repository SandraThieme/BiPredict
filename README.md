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
If you want to get only predictions for bicliques above a munimum c/p-biclique size,
you can change the values of chemical_border and protein_border (otherweise c=2, p=2 are used).
The input file should be an edge list with two columns, one for compounds and one for proteis,
with columnnames "chemical" and "ProteinID" at the head of the input file (see example 'ecoli_edgelist.txt').

## Parameters
```
EDGE_LIST: Input edge list with two columns named: "chemical" and "ProteinID"
chemical_border: Least number of compounds
protein_border: Least number of proteins
CORES: Number of cores available for parallel computing (please check how many cores you have available, or set to 1, do not use all available cores)
```

## Usage 
```
source('bipredict.R')
edge_list = read.table('ecoli_edgelist.txt',header = T)
predicted_interactions_52 = bipredict(edge_list,chemical_border = 5,protein_border = 2,CORES = 1)
write.table(predicted_interactions_52, file = 'my_predicted_interactions_52.txt', row.names = F)
```

Please find descriptions and examples also in the file bipredict.R You need both .R files in the same directory to run the predictions.


# Authors
Sandra Thieme
PhD Student
AG Bioinformatics

Max Planck Institute for Molecular Plant Physiology
Wissenschaftspark Golm, Am Mühlenberg 1, 14476 Potsdam-Golm 
https://www.mpimp-golm.mpg.de/7955/2walther
