#!/bin/bash

#Выравнивание современных людей
cat *.fasta > human.fasta
mafft --auto human.fasta > human_aligned.fasta

# Дерево с бутстрэпом
iqtree -s human_aligned.fasta -B 1000

# Добавление внешней группы - бонобо
wget "https://www.ncbi.nlm.nih.gov/kis/download/sequence/?db=nucleotide&id=NC_001644.1"
cat human_aligned.fasta *NC_001644.1* > outgroup.fasta
mafft --auto outgroup.fasta > outgroup_aligned.fasta

iqtree -s outgroup_aligned.fasta -B 1000

# Загружаем дерево в MEGA в форме .nwk и укореняем дерево по внешней группе
# Добавляем выравнивание outgroup_aligned.fasta для вычисления pairwise distance и возраста митохондриальной Евы

# Добавляем данные денисовцев и неандертальцев
cat outgroup_aligned.fasta ./Denisova/*.fasta ./Neanderthal/*.fasta > combined.fasta
mafft --auto combined.fasta > combined_aligned.fasta
iqtree -s combined_aligned.fasta -B 1000 

# Снова загружаем в MEGA и укореняем дерево

