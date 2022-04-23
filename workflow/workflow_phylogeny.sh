#!/bin/bash

# As we discussed yesterday, I have attached the sequences of Aedes Japonicus.
# We obtained using Sanger sequencing and PCR amplification of cytochrome oxidase I. 
# Please we need a plylogenetic tree and a comparison with Aedes canadensis
# https://www.ncbi.nlm.nih.gov/nuccore/?term=Aedes%20canadensis

conda activate --stack /users/PAS0471/jelmer/.conda/envs/mafft-env
conda activate --stack /users/PAS0471/jelmer/miniconda3/envs/fasttree-env

## Download Ae. canadensis sequences
conda activate /users/PAS0471/jelmer/miniconda3/envs/blast-env
esearch -db nuccore -query '"Ochlerotatus canadensis"[Organism] AND COI' |
    efetch -format fasta > data/ref/AeCanadensis_COI_nuccore.fa

## Add Iowa sequencxe
query="GQ254798.1"; esearch -db nuccore -query $query | efetch -format fasta > data/ref/"$query".fa

## Get outgroups from here https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7400179/
esearch -db nuccore -query "JX629050.1" | efetch -format fasta > data/ref/Torrenticola.fa
esearch -db nuccore -query "KX421869.1" | efetch -format fasta > data/ref/Lebertia.fa

## Align the sequences
### Wooster + Ferdinand's japonicus + 1 canadensis
#cat data/ref/canadensis_HQ944971.1.fa data/ref/ferdinand_refs.fa data/wooster.fa > results/unalign.fa
#mafft --reorder --auto --adjustdirection --leavegappyregion results/unalign.fa > results/align.fa 

### Wooster + Ferdinand's japonicus + 1 canadensis + 2 far outgroups
cat data/ref/GQ254798.1.fa data/ref/canadensis_HQ944971.1.fa data/ref/Torrenticola.fa data/ref/Lebertia.fa data/ref/ferdinand_refs.fa data/wooster.fa > results/unalign_woutgroups.fa
mafft --reorder --auto --adjustdirection --leavegappyregion results/unalign_woutgroups.fa > results/align_woutgroups.fa 

### Wooster + Ferdinand's japonicus + all canadensis + 2 far outgroups
#cat data/ref/AeCanadensis_COI_nuccore.fa data/ref/Torrenticola.fa data/ref/Lebertia.fa data/ref/ferdinand_refs.fa data/wooster.fa > results/unalign_woutgroups_allcanad.fa
#mafft --reorder --auto --adjustdirection --leavegappyregion results/unalign_woutgroups_allcanad.fa > results/align_woutgroups_allcanad.fa

## Make a tree with IQ-Tree
sbatch -t5:00:00 mcic-scripts/trees/iqtree.sh -i results/align_woutgroups.fa -p results/iqtree/woutgroups/COI -a "--msub mitochondrial -o KX421869.1" -b 1000 
Rscript scripts/plot_tree.R results/iqtree/woutgroups/COI.treefile results/iqtree/COI.png KX421869.1

## Make a tree with Fasttree
cat results/align_woutgroups.fa | FastTree -gamma -nt -gtr -out results/fasttree/woutgroups.tree
conda activate --stack base
Rscript scripts/plot_tree.R results/fasttree/woutgroups.tree results/fasttree/woutgroups.png KX421869.1



# ------------
# Include other japonicus sequences?
conda activate /users/PAS0471/jelmer/miniconda3/envs/blast-env
esearch -db nuccore -query '"Aedes japonicus"[Organism] AND COI"' -country "USA" |
    esummary | xtract -pattern DocumentSummary -element AccessionVersion SubName
# MN509209.1      DNAS-42A-5ZBI|female|adult|USA: New York, Cold Spring Harbor, Cold Spring Harbor Laboratory DNALC|Specimen obtained from behind DNALC from mosquito trap|40.87 N 73.45 W|13-Aug-2019|Austin Sarker-Youn, Jeffry Petracca, Cornel Ghiban, Sharon Pepenella, Cristina Marco, Bruce Nash, Dave Micklos|Jeffry Petracca
# MN451704.1      DNAS-40C-5XJI|female|adult|USA: New York, Cold Spring Harbor, Cold Spring Harbor Laboratory DNALC|behind DNALC from mosquito trap|40.87 N 73.45 W|17-Jul-2019|Ellie Bartolino, Jeffry Petracca, Cornel Ghiban, Sharon Pepenella, Cristina Marco, Bruce Nash, Dave Micklos|Jeffry Petracca
# MN058612.1      Everett1|larvae|USA: Everett, Washington
# MK372912.1      WaC|USA: Washington County, AL
# MK372911.1      MaC|USA: Marshall County, AL
# MK372910.1      FaC|USA: Fayette County, AL
# GQ254801.1      VT_397|USA|41.5832 N 93.6839 W|25-Jul-2008|Brendan Dunphy
# GQ254800.1      VT_624|USA|41.5430 N 90.5343 W|08-Aug-2008|Brendan Dunphy
# GQ254799.1      VT_950|USA|42.5142 N 90.7295 W|22-Aug-2008|Brendan Dunphy => Iowa
# GQ254798.1      VT_1079|USA|43.2983 N 91.8004 W|22-Aug-2008|Brad Tucker|Brendan  => Iowa
# GQ254797.1      MPINHS-G|USA|40.0843 N 88.2175 W|17-Jul-2007|Richard Lampman
# GQ254796.1      MPINHS-A|USA|38.7407 N 88.0938 W|07-May-2007|Richard Lampman
# GQ254795.1      MPIHNS-C|USA|42.0184 N 87.9981 W|24-Jul-2007|Richard Lampman => Illinois
# GQ254794.1      MPINHS-D|USA|40.0973 N 88.2265 W|31-Jul-2007|Richard Lampman => Illinois
# GQ254793.1      MPINHS-E|USA|38.7407 N 88.0938 W|07-May-2007|Richard Lampman => Illinois
# JX259646.1      NEONTculicid2051|USA: Massachusetts, Worcester County, Harvard Forest|42.544 N 72.176 W|30-Jun-2010|Angelo, S & Grace, C|M. J. Weissmann
# JX259645.1      NEONTculicid1551|USA: Massachusetts, Worcester County, Harvard Forest|42.531 N 72.189 W|14-Jul-2010|Angelo, S & Grace, C|R. F. Darsie Jr.
# JX259644.1      NEONTculicid1535|USA: Massachusetts, Worcester County, Harvard Forest|42.544 N 72.176 W|29-Jun-2010|Angelo, S & Grace, C|R. F. Darsie Jr.
# JX259643.1      NEONTculicid2055|USA: Massachusetts, Worcester County, Harvard Forest|42.544 N 72.176 W|30-Jun-2010|Angelo, S & Grace, C|M. J. Weissmann
# JX259642.1      NEONTculicid2054|USA: Massachusetts, Worcester County, Harvard Forest|42.544 N 72.176 W|30-Jun-2010|Angelo, S & Grace, C|M. J. Weissmann
# JX259641.1      NEONTculicid2053|USA: Massachusetts, Worcester County, Harvard Forest|42.544 N 72.176 W|30-Jun-2010|Angelo, S & Grace, C|M. J. Weissmann
# JX259640.1      NEONTculicid2052|USA: Massachusetts, Worcester County, Harvard Forest|42.544 N 72.176 W|30-Jun-2010|Angelo, S & Grace, C|M. J. Weissmann
# JX259639.1      NEONTculicid2056|USA: Massachusetts, Worcester County, Harvard Forest|42.544 N 72.176 W|30-Jun-2010|Angelo, S & Grace, C|M. J. Weissmann
# HQ978778.1      BIOUG<CAN>:TDWG-0240|USA: Massachusetts|41.5299 N 70.652 W|27-Sep-2010|BIObus: A.Borisenko, R.Labbee, V.Levesque-Beaudin, M.Milton, J.Smith, S.Ratnasingham
# HQ978777.1      BIOUG<CAN>:TDWG-0239|USA: Massachusetts|41.5299 N 70.652 W|27-Sep-2010|BIObus: A.Borisenko, R.Labbee, V.Levesque-Beaudin, M.Milton, J.Smith, S.Ratnasingham