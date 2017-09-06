# Introduction

Welcome to the Temporal Network Centralities (TNC) package. It has been
developed to calculate temporal centrality values for nodes of dynamic networks represented as a sequence of static graphs.
At the moment there are three measures available that use a snapshot based representation of the network:

* temporal betweenness centrality (TBC)

* temporal closeness centrality (TCC)

* temporal degree centrality (TDC)

For a formal definition of the measures see e.g. the paper *Temporal node centrality in complex networks* (Kim \& Anderson, Physical Review E, 2012).

While TDC is simply the average of the degree centrality of every node at every snapshot, calculating TBC and TBC is more tricky. These measures are based on **temporal shortest paths** that have to be calculated not only once but for all possible trunkated snapshot sequences of the original snapshot sequence. Hence, the computation can easily become time and memory demanding with a growing number of nodes and snapshots. To reduce the computational burden we have implemented the **Reversed Evolution Network (REN)** algorithm (see Hanke \& Foraita, *Clone temporal centrality measures for incomplete sequences of graph snapshots*, BMC Bioinformatics, 2017). REN's computational effort is linear in the number of snapshots. However, the calculation of TCC and TBC is still quadratic and cubic in the number of nodes, respectively. We have therefore implemented the possibility to parallelize the calculation of TBC and TCC by searching only for temporal shortest paths ending at a user specified set of nodes. 

Note, the TNC package accepts both adjacency lists and adjacency matrices as input. Matrices, however, exponentially increase the required memory even if the dynamic temporal network is sparse. Hence, we suggest to use adjacency lists to calculate temporal centrality values. 

# Loading
```{r}
library(TNC)
```

# Usage
First we create a toy list of adjacency matrices:
```{r}
A1 <- matrix(c(0,1,0,0,0,0,
               1,0,1,0,0,0,
               0,1,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,0,0,0), ncol=6)
A2 <- matrix(c(0,0,0,0,0,0,
               0,0,1,0,0,0,
               0,1,0,1,1,0,
               0,0,1,0,0,0,
               0,0,1,0,0,0,
               0,0,0,0,0,0), ncol=6)
A3 <- matrix(c(0,0,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,0,0,0), ncol=6)
A4 <- matrix(c(0,1,0,0,0,0,
               1,0,0,1,0,0,
               0,0,0,0,0,0,
               0,1,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,0,0,0), ncol=6)

As <- list(A1,A2,A3,A4)
```
Next, we calculate TBC, TCC and TDC:
```{r}
tbc(As, type="M")
tcc(As, type="M")
tdc(As, type="M")
```
Since the temporal centrality of a node can vary along the snapshot sequence, the centrality can be monitored over the time span setting `centrality_evolution = TRUE`:
```{r}
tbc(As, type="M", centrality_evolution = TRUE)
```
The output is a list where the first element is a vector of the temporal centrality values for each node and the second elements is a p x T matrix containing in each column the temporal centrality values for each node measured at start snapshot 1,2,...,T of the subsequences (this is the temporal centrality at each snapshot). 
