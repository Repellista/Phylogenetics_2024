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

# Загружаем дерево в MEGA в форме .nwk и укореняем дерево по внешней группе. При необходимости нажимаем Toggle the tree size
# Добавляем выравнивание outgroup_aligned.fasta для вычисления pairwise distance и возраста митохондриальной Евы

# Добавляем данные денисовцев и неандертальцев
cat outgroup_aligned.fasta ./Denisova/*.fasta ./Neanderthal/*.fasta > combined.fasta
mafft --auto combined.fasta > combined_aligned.fasta
iqtree -s combined_aligned.fasta -B 1000 

# Снова загружаем в MEGA и укореняем дерево
# Открываем TimeTree и загружаем выравнивание и дерево (при необходимости меняя названия)
# Выделяем внешнюю группу
# Находим в графическом интерфейсе узел, соответствующий предку современных людей (от которого отделяется центральноафриканская линия)
# Нажимаем на значок калибровки и указываем Fixed Date (174054). Получаем калиброванное дерево.
