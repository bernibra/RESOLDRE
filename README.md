##muritz

**RESOLDRE** is an R package for generating “correlation-informed” null models, which combine the classic concept of null models and tools from joint statistical modeling in community ecology. Such models can be used to assess whether the information encoded within any given correlation matrix is predictive for explaining structural patterns observed within an incidence matrix.


## How to cite?
To use this software, please make sure you cite Bramon Mora et. al. (*Unmasking common structural patterns in ecological data. Journal*, year).

## Installation instructions

Installation should be relatively painless:

you must clone the repository

		git clone git@github.com:bernibra/RESOLDRE.git

and build and install the package using R.


## Data accessability

The data used in Bramon Mora et. al. (year) to present the software can be find in the './data' directory. In particular, we used two different datasets:

1. The first example used in the paper is an application of the software to food webs in order to study how well species' evolutionary relationships can explain observed patterns of predator-pray interactions. The data used for this can be found in './data/foodwebs' and describes 10 empirical food webs from small streams of the Taieri River in New Zealand comprising fish, macroinvertebrates and algae (Townsend et al., 1998). Every directory in './data/foodwebs' describes a different food web (in total, there are 10 food webs). For every food web, we provide three files 'interactions.txt', 'species.txt' and 'output_tree.new' describing the list of of interactions composing the food webs, the list of species involved in these interactions and the phylogenetic tree characterizing the evolutionary history of such species, respectively. Notice that the code used for the different species in 'interactions.txt' is an integer number that represents the position of the species in 'species.txt'. Although the data is available online, make sure you cite Townsend et al. (1998) if you want to use it.

1. The second example used in the paper is an application of the to software species assemblages in order to study how well possible spatial autocorrelations or area similarity between sample sites as well as island richness and species range similarity can explain the structural patterns observed in these communities. The data used for this can be found in './data/biogeography' and describes the distribution of 366 species of vascular plants across 80 islands from the San Juan archipelago (Marx et al., 2015a; Marx et al., 2015b). In particular, we provide two files 'DRYAD1\_ComMatrix.csv' and 'DRYAD2\_SJtraits.csv' describing the incidence matrix for this species assemblage and the list of species traits, respectively.  Although the data is available online, plase make sure you cite Marx et al. (2015a; 2015b) if you want to use it.
