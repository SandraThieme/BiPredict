# BiPredict
Biclique based compound-protein-interaction (CPI) prediction


Bicliques consist of two types of nodes and edges connecting each node of different types (a). Here, blue circles represent compounds, c, and red squares represent proteins, p, while edges represent interactions. The biclique expansion starts with an existing biclique (b, yellow inner circle), here consisting of three compounds (c = 3) and two proteins (p = 2). Next, all compounds and proteins that are directly connected to any member of the biclique are identified. They represent interaction candidates (b, lightblue outer circle). Other compounds and proteins of the network which are not directly connected to any member of an existing biclique are not considered (b, grey squares and circles). Finally, interactions are predicted if interaction candidates lack only one edge to be a member of an existing biclique (b, green dashed lines). 
![methods_skizze](https://user-images.githubusercontent.com/82212543/126770821-6b673a8d-7bdc-4036-b748-0adf11509ded.png)

Please find descrition and examples in the file bipredict.R You need both .R files in the same directory to run the predictions. 
You also need to install the following dependencies (R packages):
biclique,
igraph,
parallel,
data.table

# Authors
Sandra Thieme
PhD Student
AG Bioinformatics

Max Planck Institute for Molecular Plant Physiology
Wissenschaftspark Golm, Am Mühlenberg 1, 14476 Potsdam-Golm 
https://www.mpimp-golm.mpg.de/7955/2walther
