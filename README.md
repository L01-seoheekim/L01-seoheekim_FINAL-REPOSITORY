# final-project-L01-seoheekim
Below is the code that is required to regenerate the analysis of the SLC gene family specifically on the SLC36 Nematostella vectensis query sequence of XP_032223257.1. The code will be divided in four sections that in the order in which they should be inputed.

# SLC Identification: Sequencing and Alignment

```
ncbi-acc-download -F fasta -m protein XP_032223257.1
#Download of query sequence from NCBI

blastp -db allprotein.fas -query XP_032223257.1.fa -outfmt 0 -max_hsps 1 -out AAT.blastp.typical.out
#Performing of blast search using query protein

less AAT.blastp.typical.out

blastp -db allprotein.fas -query XP_032223257.1.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle"  -max_hsps 1 -out AAT.blastp.detail.out
#More detailed blast analysis

less AAT.blastp.detail.out

awk '{if ($6<0.00000000000001)print $1 }' AAT.blastp.detail.out > AAT.blastp.detail.filtered.out
#Filtering of putative homologs, setting e-value less than 1e -14

wc -l AAT.blastp.detail.filtered.out 

seqkit grep --pattern-file AAT.blastp.detail.filtered.out allprotein.fas > AAT.blastp.detail.filtered.fas
#Change file format from .out to .fas

muscle -in AAT.blastp.detail.filtered.fas -out AAT.blastp.detail.filtered.aligned.fas
#Alignment of sequence 

t_coffee -other_pg seq_reformat -in AAT.blastp.detail.filtered.aligned.fas -output sim
#To provide statistics of alignment

alv -k AAT.blastp.detail.filtered.aligned.fas | less -r
#To view alginment

alv -kli --majority AAT.blastp.detail.filtered.aligned.fas | less -RS
#To view alignment

t_coffee -other_pg seq_reformat -in AAT.blastp.detail.filtered.aligned.fas -action +rm_gap 50 -out AAT.blastp.detail.filtered.aligned.r50.fas
#Removal of highly gapped regions, specifically gaps greater than 50% gapped residues

alv -kli --majority AAT.blastp.detail.filtered.aligned.r50.fas | less -RS
#To view alignment
```
 
 # Formation of Phylogenetic Tree
 
 ```
 sed "s/ /_/g" AAT.blastp.detail.filtered.aligned.fas > AAT.blastp.detail.filtered.aligned_.fas
 # Replacing the indication of spaces to underscores to prevent lossing space annotations
 
 iqtree -s AAT.blastp.detail.filtered.aligned_.fas -nt 2
 # To get maximum likehood tree estimate
 
 gotree reroot midpoint -i AAT.blastp.detail.filtered.aligned_.fas.treefile -o AAT.blastp.detail.filtered.aligned_.fas.midpoint.treefile
 # Rerooting of tree to midpoint

nw_display -s  maguk.blastp.detail.filtered.aligned_.fas.midpoint.treefile -w 1000 -b 'opacity:0' >  maguk.blastp.detail.filtered.aligned_.fas.midpoint.treefile.svg
# To view the created tree
 ```
 
 # Reconciliation of Rooted Tree
 
 ```
 java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -s speciesTreeBilateriaCnidaria.tre -g AAT.blastp.detail.filtered.aligned_.fas.midpoint.treefile --reconcile --speciestag prefix --savepng --events
 # Reconciliation of species tree with rerooted midpoint tree
 
 python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g AAT.blastp.detail.filtered.aligned_.fas.midpoint.treefile.reconciled --include.species
 # Convertion of format to .xml file
 
 thirdkind -f AAT.blastp.detail.filtered.aligned_.fas.midpoint.treefile.reconciled.xml -o  AAT.blastp.detail.filtered.aligned_.fas.midpoint.treefile.genes.tre.reconciled.svg
 # Reformtating of reconciled tree from .xml to .svg file
 
 java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -s speciesTreeBilateriaCnidaria.tre -g AAT.blastp.detail.filtered.aligned_.fas.midpoint.treefile --root --speciestag prefix --savepng --events
 # Minimizes duplications and deletions in reconciled tree
 
iqtree -s AAT.blastp.detail.filtered.aligned_.fas -bb 1000 -nt 2 -m LG+F+R5 -t AAT.blastp.detail.filtered.aligned_.fas.midpoint.treefile -pre AAT.genes.ufboot
 # The running of iqtree to preform bootstrap to obtain the most optimal tree with highest score

gotree reroot midpoint -i AAT.blastp.detail.filtered.aligned_.fas.treefile -o AAT.genes.midpoint.ufboot
 # Rerooting of tree to correspond with caluclated bootstrap
 
 java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -s speciesTreeBilateriaCnidaria.tre -g AAT.blastp.detail.filtered.aligned_.fas.midpoint.treefile.rooting.0 --root --speciestag prefix --savepng --events
 # Minimizes duplications and deletions in new reconciled tree, rerooting it
 ```
 
 # Domain Identification in SLC Protein
 
 ```
 iprscan5   --email seohee.kim.1@stonybrook.edu  --multifasta --useSeqId --sequence   AAT.blastp.detail.filtered.fas
 # Sends query protein to interproscan server
 
 cat ~/labs/lab8-L01-seoheekim/myfamily/*.tsv.tsv > ~/labs/lab8-L01-seoheekim/AAT.domains.all.tsv
 # Concentrates the identified domains for each gene into one .tsv file
 
 grep Pfam ~/labs/lab8-L01-seoheekim/AAT.domains.all.tsv >  ~/labs/lab8-L01-seoheekim/AAT.domains.pfam.tsv
 # Filterning of domains not identified in Pfam database
 
 awk 'BEGIN{FS="\t"} {print $1"\t"$3"\t"$7"@"$8"@"$5}' ~/labs/lab8-L01-seoheekim/AAT.domains.pfam.tsv | datamash -sW --group=1,2 collapse 3 | sed 's/,/\t/g' | sed 's/@/,/g' > ~/labs/lab8-L01-seoheekim/AAT.domains.pfam.evol.tsv
 # Rearrangement of interproscan output
```
