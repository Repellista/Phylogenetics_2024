# Introduction
This research is aimed at clarifying the phylogenetic position of different Testudines, namely Chelonidae, Dermochelyidae, Testudinidae, Emydidae and Trioychidae families. Chelonidae and Dermochelyidae belong to sea turtles from a common superfamily Chelonioidea, whereas Testudinidae include terrestial tortoises. Emydidae is represented by asemiaquatic Trachemys scripta elegans, and Trionyx triunguis is a freshwater turtle from Trioychidae family. The research is aimed as clarifying the order of divergence of sea, terrestrial and freshwater turtles. Study[1] claimed that Emydidae were more closely related to Testudinidae. Study[2] was telling that Tryonichia stand distinct to other turtles.

# Methods
Samples for 11 species of turtles were accessed from Genbank. Table with accession numbers was taken from study [2] and some other species were added for better family resolvance. Alligator missipiensis was chosen as outgroup. 
```
import Bio
from Bio import Entrez
from Bio import SeqIO
Entrez.email = "..."
import pandas as pd
table=pd.read_csv("/content/tortoise_data_numbers.csv", sep=";")
lst12S = table[table['Gene']=='12S']
an_12S = lst12S['Accession number'].str.split('[,;]', expand=False).explode().tolist()
handle1 = Entrez.efetch(db="nucleotide", id=an_12S, rettype="fasta", retmode="text")
filename='12S.fasta'
with open(filename, 'w') as file:
  file.write(handle1.read())
```
In the same manner, data for other genes and Alligator missipiensis was retrived.

Alignment was done with MAFFT for each gene marker. 
```
mafft --auto turtles.fasta > turtles_aligned.fasta
```

Muscle was also used for alignment, but it didn't result in any difference.
```
muscle -in turtles.fasta -out turtles_muscle_aligned.fasta
```

Several trees were built with different approaches and visualized with ITOL: neighbor-joining with phylip, maximum likelihood with IQTree and Bayesian with MrBayes.

The commands used for NJ approach:
```
seqret -sequence turtles_aligned.fasta -outseq turtles_aligned.phy -osformat2 phylip
dnadist
neighbor
```
F84 was used as metric and 0.5 as alpha.

Maximum Likelihood with IQTree:
```
iqtree -s turtles_aligned.fasta -m GTR+G -bb 1000 -alrt 1000
```
MrBayes:
```
seqret -sequence turtles_aligned.fasta -outseq turtles_aligned.nex -osformat2 nexus

mb
execute supermatrix.nex
mcmc ngen=1000000 samplefreq=100 printfreq=1000 nchains=4
sump burnin=250
sumt burnin=250
```
# Results: 
NJ, ML and Bayesian methods all resulted in trees with identical topology and high bootstrap support = 100 %.

![00eWQegcnEePg1kbYKzrNA](https://github.com/user-attachments/assets/8027a27a-5294-4f02-96cd-d638cdb8a69a)

# Discussion 
The topological order derived alignes study [3]. The first to diverge were the members of Trioychidae family, softshell turtles whose inhabitat is semiaquatic waters. 
Then, the lineage splits into two parallel lineages, one being representatives of Chelonioidea superfamily, namely sea turtles, and another node which includes both 
terrestrial tortoises and Emydidae family, who belong to freshwater. We can conclude this way that modern terrestrial tortoises developed their lifestyles later 
than water turtles. Also, freshwater turtles are closer to terrestrial than sea ones.

# Literature
1. "Molecular phylogenetics of some endangered turtles reveals new close genetic relationships" ( https://doi.org/10.26577/ijbch.2023.v16.i1.02)
2. "Phylogenomic analyses of 539 highly informative loci dates a fully resolved time tree for the major clades of living turtles (Testudines)" (https://doi.org/10.1016/j.ympev.2017.07.006)
3. A global phylogeny of turtles reveals a burst of climate-associated diversification on continental margins  (doi: 10.1073/pnas.2012215118)
