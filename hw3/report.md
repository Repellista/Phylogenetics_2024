Introduction
This research is aimed at clarifying the phylogenetic position of different Testudines, namely Chelonidae, Dermochelyidae, Testudinidae, Emydidae and Trioychidae families. 
Chelonidae and Dermochelyidae belong to sea turtles from a common superfamily Chelonioidea, whereas Testudinidae include terrestial tortoises. Emydidae is represented by asemiaquatic Trachemys scripta 
elegans, and Trionyx triunguis is a freshwater turtle from Trioychidae family. 
The assumptions were that Testudinidae, being the only terrestrial turtles here, form a separate clade to all others. Study[1] claimed that Emydidae were more closely related to Testudinidae, and Chelonioidea formed sequential sister groups to them, with Dermochelyidae forming a monophyletic group. Study[2] opposed it by telling that Tryonichia stand distinct to other turtles and that Dermochelyidae is a sister taxon to Chelonidae.

Methods
Samples for 9 species of turtles were accessed from Genbank. Alligator missipiensis was chosen as outgroup. Four gene markers were chosen: 12S, 16S, COI and RAG-1.
As there were several sequences for the same gene marker and species, consensus sequences were made for them and aligned with MAFFT for each gene marker. Then, supermatrix of concatenated genes was built with AMAS.
Both separate gene markes and concatened supermatrix were used as data, as well as all sequences ungrouped.
Several trees were built with different approaches and visualized with ITOL: neighbor-joining with phylip, maximum likelihood with IQTree and Bayesian with MrBayes.

Results
There were significant differences in inferred phylogenetic relationships between species when using different markers. 
For ML approach, the most reliable bootstrap support values were found for 12S gene, others had weak support. 

12S
![dTMTDWC29BssgV9ylyyuzA](https://github.com/user-attachments/assets/f2e3979a-381e-4404-a041-81aa913cb344)

16S
![pJcKuURWA4Vy2oI1YbB6Dw](https://github.com/user-attachments/assets/5e9e49c4-5abe-4f96-8fec-5d8558cd89a4)

COI
![uRF1MEVyQ4lY-XgETepkSQ](https://github.com/user-attachments/assets/4a150c9a-1807-4178-8502-db38f2ec2d79)

RAG-1
![IfSzGHLR-nu1mCVEjkzGnA](https://github.com/user-attachments/assets/fadbc4cc-57a2-40db-80cc-d8d750cbd6a4)


As for NJ trees, the 12S tree is consistant with ML one, whereas others also have weak support and negative branch length values.

![kk5hvX4rrhbv94tft1W99A](https://github.com/user-attachments/assets/09f8f90d-5b31-4760-915f-3c2fdbd30f48)

Concatenated NJ tree: 

![NJ_outgroup_consensus_001](https://github.com/user-attachments/assets/9fd0c802-9695-4ae4-a1ef-f54dabbb1c76)

Concatenated ML tree:

![ML_consensus](https://github.com/user-attachments/assets/77e05157-f381-4a81-a97e-9591b63cba95)

Concatenated Bayes tree:

![Bayes_consensus](https://github.com/user-attachments/assets/c9b6a693-d377-478b-b8b3-6b7cdc3bad6a)

The tree with all sequences was also made to check some of the reasons for weak support:

![all_sequences](https://github.com/user-attachments/assets/cb58b5f1-8bbc-44f6-a914-1bf4c49fa6d4)


Discussion
The derived trees do not share consistency with each other. Judging by most trees, Testudines do stand furthest from other turtles. That has sense, since they have different lifestyle than other turtles that live in water. However, the order of separating is reversed for NJ and ML method, where NJ suggests Testudines were the first to separate from the outgroup, and ML suggests they were the last. 
Dermochelys coriaces mostly stands close to Chelonidae representatives, but clearly forms a separate group, proving that they are from different families but the same superfamily.
As for Emydidae group, in most trees, it is close to Testudines and suggest that it does not form any group distinct from all others. However, Bayesian tree based on supermatrix data, indeed suggests that Trachemys scripta diverges the first, however the branching order with Tryonix triunguis is not clear due to polytomy. 
Also, the findings oppose the claim in study [2] that Dermochelys coriaces is more closely related to Caretta caretta than Chelonia mydas when based on RAG-1 and COI genes. 
Overall, supermatrix-based ML tree repeats the scheme of branching order described in study [1].  

Literature
1. "Molecular phylogenetics of some endangered turtles reveals new close genetic relationships" ( https://doi.org/10.26577/ijbch.2023.v16.i1.02)
2. "Phylogenomic analyses of 539 highly informative loci dates a fully resolved time tree for the major clades of living turtles (Testudines)" (https://doi.org/10.1016/j.ympev.2017.07.006)
