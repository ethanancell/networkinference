---
title: "Post-selection inference with a single realization of a network"
---

## Post-selection inference with a single realization of a network

The `networkinference` R package implements the methods from ``Post-selection inference with a single realization of a network'' by Ethan Ancell, Daniela Witten, and Daniel Kessler.

Suppose that you observe a network of $n$ nodes, encoded in the adjacency matrix $A \in \mathbb{R}^{n \times n}$ where $A_{ij}$ encodes the edge from node $i$ to node $j$. If the network is unweighted, then $A_{ij} = 1$ if there is an edge from node $i$ to node $j$, and $A_{ij} = 0$ otherwise. If the network is weighted, then $A_{ij}$ is the weight from the edge pointing from node $i$ to node $j$. If the network is undirected, then $A_{ij} = A_{ji}$ for all $i$ and $j$.

Our setting is that of *random networks*: in particular, we assume one of the following holds:

* $A_{ij} \overset{\text{ind.}}{\sim} N(M_{ij}, \sigma^2)$, where $\sigma^2$ is known
* $A_{ij} \overset{\text{ind.}}{\sim} \text{Poisson}(M_{ij})$
* $A_{ij} \overset{\text{ind.}}{\sim} \text{Bernoulli}(M_{ij})$

This package and the associated paper considers the setting where the user wishes to 

(1) Estimate latent communities in the network based upon the adjacency matrix.
(2) Conduct inference for a connectivity parameter that is a function of those estimated communities. In particular, we consider a connectivity parameter that is a linear combination of the average connectivity within and between the estimated communities.

For step (2) to be valid, it *must* take into account that the selected parameter is a function of the data, and so the coverage of confidence intervals must hold conditional on the selection of the communities.

To accomplish this, we "split" an adjacency matrix $A \in \mathbb{R}^{n \times n}$ into two adjacency matrices $A^{(\text{tr})}, A^{(\text{te})} \in \mathbb{R}^{n \times n}$ using the `networkinference::split_network()` function. Then, the user may use any method of their choice to estimate communities $\hat{Z} = \hat{Z}(A^{(\text{tr})}) \in \{0,1\}^{n \times K}$ where $\hat{Z}_{ik} = 1$ if the $i$th node belongs to the $k$th estimated community, and $\hat{Z}_{ik} = 0$ otherwise. Critically, $\hat{Z}$ must only be a function of $A^{(\text{tr})}$.

Finally, we conduct inference for the selected parameter using the `networkinference::infer_network()` function.

## Examples

For an end-to-end example running through all the features in this R package, please read [this vignette](articles/network-inference-vignette.html).

## Installation

To install the `networkinference` package, run `devtools::install_github("ethanancell/networkinference")` in your R console.

## Figures in paper

[Here is a link to the Github repository](https://github.com/ethanancell/networkinferencepaper) that contains all of the scripts used to generate the figures in the paper.

