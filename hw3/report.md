# Introduction
This research is aimed at clarifying the phylogenetic position of different Testudines, namely Chelonidae, Dermochelyidae, Testudinidae, Emydidae and Trioychidae families. 
Chelonidae and Dermochelyidae belong to sea turtles from a common superfamily Chelonioidea, whereas Testudinidae include terrestial tortoises. Emydidae is represented by asemiaquatic Trachemys scripta 
elegans, and Trionyx triunguis is a freshwater turtle from Trioychidae family. 
The assumptions were that Testudinidae, being the only terrestrial turtles here, form a separate clade to all others. Study[1] claimed that Emydidae were more closely related to Testudinidae, and Chelonioidea formed sequential sister groups to them, with Dermochelyidae forming a monophyletic group. Study[2] opposed it by telling that Tryonichia stand distinct to other turtles and that Dermochelyidae is a sister taxon to Chelonidae.

# Methods
Samples for 9 species of turtles were accessed from Genbank. Table with accession numbers was taken from study [2]. Alligator missipiensis was chosen as outgroup. Four gene markers were chosen: 12S, 16S, COI and RAG-1.
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

As there were several sequences for the same gene marker and species, consensus sequences were also made for them.

```
import os
from Bio import SeqIO
from Bio.Align.Applications import MafftCommandline
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
import subprocess

def extract_species(header):
    parts = header.split()
    species = parts[1] + " " + parts[2]  
    return species

def group_sequences_by_species(input_fasta):
    """Group sequences by species based on their FASTA headers."""
    species_sequences = defaultdict(list)
    for record in SeqIO.parse(input_fasta, "fasta"):
        species = extract_species(record.description)
        species_sequences[species].append(record)
    return species_sequences

def write_species_fasta(species, records, output_dir):
    """Write grouped sequences to individual FASTA files per species."""
    output_fasta = os.path.join(output_dir, f"{species.replace(' ', '_')}.fasta")
    with open(output_fasta, "w") as out_f:
        SeqIO.write(records, out_f, "fasta")
    return output_fasta

def run_mafft(input_fasta, output_fasta):
    """Align sequences using MAFFT."""
    mafft_cline = MafftCommandline(input=input_fasta)
    stdout, stderr = mafft_cline()
    with open(output_fasta, "w") as out_f:
        out_f.write(stdout)
    return output_fasta

def generate_consensus(aligned_fasta, species):
    """Generate a consensus sequence from an aligned FASTA file."""
    alignment = AlignIO.read(aligned_fasta, "fasta")
    consensus = ""
    for i in range(len(alignment[0].seq)):
        column = alignment[:, i]  # Extract each column from alignment
        most_common_base = max(set(column), key=column.count)
        consensus += most_common_base
    return SeqRecord(Seq(consensus), id=species.replace(' ', '_'), description="consensus sequence")

def main(input_fasta, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    species_sequences = group_sequences_by_species(input_fasta)
    
    consensus_records = []
    
    for species, records in species_sequences.items():
        print(f"Processing species: {species}")
        
        # Write sequences to species-specific FASTA
        species_fasta = write_species_fasta(species, records, output_dir)
        
        # Align sequences using MAFFT
        aligned_fasta = os.path.join(output_dir, f"{species.replace(' ', '_')}_aligned.fasta")
        run_mafft(species_fasta, aligned_fasta)
        consensus_record = generate_consensus(aligned_fasta, species)
        consensus_records.append(consensus_record)
    
    # Write all consensus sequences to a single FASTA file
    output_fasta = os.path.join(output_dir, "consensus_sequences.fasta")
    with open(output_fasta, "w") as out_f:
        SeqIO.write(consensus_records, out_f, "fasta")
    
    print(f"Consensus sequences saved to {output_fasta}")


input_fasta = "input_sequences.fasta" 
output_dir = "consensus_output"
main(input_fasta, output_dir)
```

and aligned with MAFFT for each gene marker. 
```
mafft --auto 12S_consensus.fasta > 12S_consensus_aligned.fasta
```

Muscle was also used for alignment, but it didn't result in any difference.
```
muscle -in 12S.fasta -out 12S_muscle_aligned.fasta
```

Then, supermatrix of concatenated genes was built with AMAS.
```
amas concat -i 12S_aligned.fasta COI_aligned.fasta RAG1_aligned.fasta -f fasta -d dna -t supermatrix.fasta -p partitions.txt
```

Both separate gene markes and concatened supermatrix were used as data, as well as all sequences ungrouped.
Several trees were built with different approaches and visualized with ITOL: neighbor-joining with phylip, maximum likelihood with IQTree and Bayesian with MrBayes.

The commands used for NJ approach:
```
seqret -sequence supermatrix.fasta -outseq supermatrix.phy -osformat2 phylip
dnadist
neighbor
```
F84 was used as metric and 0.5 as alpha.

Maximum Likelihood with IQTree:
```
iqtree -s supermatrix.fasta -p partitions.txt -m GTR+G -bb 1000 -alrt 1000
```
MrBayes:
```
seqret -sequence supermatrix.fasta -outseq supermatrix.nex -osformat2 nexus
# In the .nex file:

begin mrbayes;
    charset 12S = 1-500;     
    charset 16S = 501-1000;  
    charset COI = 1001-1500; 
    charset RAG1 = 1001-1500; 


    partition by_gene = 4: 12S, 16S, COI, RAG1;  [Combining the partitions]
    set partition=by_gene;
    
    lset applyto=(1) nst=6 rates=gamma;  
    lset applyto=(2) nst=6 rates=gamma;  
    lset applyto=(3) nst=6 rates=gamma;  
    
    prset applyto=(all) ratepr=variable;
    
    mcmc ngen=1000000 samplefreq=100 printfreq=1000 nchains=4;
    sump burnin=250;  
    sumt burnin=250;
end;

mb
execute supermatrix.nex
```


# Results
There were significant differences in inferred phylogenetic relationships between species when using different markers and phylogenetic approaches. Bayesian method showed the most bootstrap support when building trees for separate genes, but the supermatrix tree has unresolved polytomies, probably to many differences between gene models.

For ML approach, the most reliable bootstrap support values were found for 12S gene, others had weak support. 

12S
![dTMTDWC29BssgV9ylyyuzA](https://github.com/user-attachments/assets/f2e3979a-381e-4404-a041-81aa913cb344)

16S
![pJcKuURWA4Vy2oI1YbB6Dw](https://github.com/user-attachments/assets/5e9e49c4-5abe-4f96-8fec-5d8558cd89a4)

COI
![uRF1MEVyQ4lY-XgETepkSQ](https://github.com/user-attachments/assets/4a150c9a-1807-4178-8502-db38f2ec2d79)

RAG-1
![IfSzGHLR-nu1mCVEjkzGnA](https://github.com/user-attachments/assets/fadbc4cc-57a2-40db-80cc-d8d750cbd6a4)

Bayesian trees for single genes:

12S
![72eiOh4XHr5MD3xHKdyhlg](https://github.com/user-attachments/assets/b8f08d54-9164-4f2c-865c-b3f57b32a826)

16S

![cSwpBfThU-AhQh40xDTrew](https://github.com/user-attachments/assets/4c668303-6779-4d7b-98c3-02333309aabf)

COI

![MJrqah-pQ8p-jfy9Qwfp3Q](https://github.com/user-attachments/assets/d79f3ae2-9749-49eb-9dc2-1064654c280c)

RAG

![Au5UqGjj3y8-ru-qmuNTxA](https://github.com/user-attachments/assets/e8804df4-19bf-4344-b70a-91f25df6f66c)


As for NJ trees, the 12S tree is consistant with ML one, whereas others also have weak support and negative branch length values.

![kk5hvX4rrhbv94tft1W99A](https://github.com/user-attachments/assets/09f8f90d-5b31-4760-915f-3c2fdbd30f48)

Concatenated NJ tree: 

![NJ_outgroup_consensus_001](https://github.com/user-attachments/assets/9fd0c802-9695-4ae4-a1ef-f54dabbb1c76)

Concatenated ML tree:

![ML_consensus](https://github.com/user-attachments/assets/77e05157-f381-4a81-a97e-9591b63cba95)

Concatenated Bayes tree:

![Bayes_consensus](https://github.com/user-attachments/assets/c9b6a693-d377-478b-b8b3-6b7cdc3bad6a)



# Discussion
The derived trees do not share consistency with each other. The order of separating for Testudines from other turtles is reversed for NJ and ML method, where NJ suggests Testudines were the first to separate from the outgroup, and ML suggests they were the last. The same is true for analysis of separate genes, where 16S indicates that Testudines were the last to branch off from water turtles and others, especially 12S, suggest the opposite.

Dermochelys coriaces mostly stands close to Chelonidae representatives, but clearly forms a separate group, proving that they are from different families but the same superfamily. Also, the findings oppose the claim in study [2] that Dermochelys coriaces is more closely related to Caretta caretta than Chelonia mydas when based on RAG-1 and COI genes. 

As for Emydidae group, in ML and NJ trees, it is close to Testudines and suggest that it does not form any group distinct from all others. However, Bayesian tree based on supermatrix data, indeed suggests that Trachemys scripta diverges the first, however the branching order with Tryonix triunguis is not clear due to polytomy. Bayesian trees for 16S, COI and especially RAG genes also suggest that Emydidae stand furthest from other turtles. It can be explained by America being its prime areal, whereas other turtles inhabit other continents or have wider geography in general. 

Overall, supermatrix-based ML tree repeats the scheme of branching order described in study [1].  

# Literature
1. "Molecular phylogenetics of some endangered turtles reveals new close genetic relationships" ( https://doi.org/10.26577/ijbch.2023.v16.i1.02)
2. "Phylogenomic analyses of 539 highly informative loci dates a fully resolved time tree for the major clades of living turtles (Testudines)" (https://doi.org/10.1016/j.ympev.2017.07.006)
