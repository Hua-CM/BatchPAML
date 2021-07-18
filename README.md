# BatchPAML

## Introduction
This is the script for identifying positive selection on a specific branch for batches of gene families using PAML branch-site model.

## requirement
~~~shell
Biopython == 1.79
paml == 4.9j
macse == v2.04
Pandas == 1.24
~~~
*ps:* Other versions should work as well

## How to use

~~~shell
python BatchPAML.py -h
~~~

### Pattern "universe"
For universal use, you need provide a file contain two columns **without header**:

| Family_name | MSA_file |
| :---: | :---: |
| Family1 | Family1_aligned.fasta |
| Family2 | Family2_aligned.fasta |

I recommend using [macse](https://bioweb.supagro.inra.fr/macse/)

### Pattern "orthofinder"
This pattern was designed for the single copy family identified in [OrthoFinder](https://github.com/davidemms/OrthoFinder).
You just need to prepare a fasta file contain all corresponding cds sequences for the protein sequences used in OrthoFinder and specify some results file from OrthoFinder.  
*ps:* The species tree constructed by OrthoFinder is OK. There is no need to use gene tree for each gene family.

### output
The result file contain two columns **without header**

| Family_name | p-value |
| :---: | :---: |
| Family1 | 1.0 |
| Family2 | 0.5921 |

## Advantages
1. convert MSA in fasta format to paml format automatically.  
2. unroot the rooted tree automatically
3. Multi thread parallel
4. Allow specify the type of codon table manually (using NCBI No.)

## Notes
1. You need mark the foreground branch in the tree manually using "#1". **Please do not insert space between species name and the marker**  
example:   
(((Human#1, chimpanzee),Fish),Fly);✓  
(((Human #1, chimpanzee),Fish),Fly);✗
   
2. Though this script search the paml bin in PATH and should be cross-platform, I recommand you specify the path of binary manually. 

3. The multiple sequence alignment (MSA) **must be in fasta format**.

4. if the MSA has potential frameshift, then this family will be skipped.

5. The process file are stored in 'BatchPAML_Results' in working directory

## Future
1. Allow MSA file in paml format directly
2. Fix frameshift problem
3. More flexible file storage path
4. Maybe a part of a  comparative genome analysis pipeline

## Contact
If you have any problem or advice pleas feel free to contact me by njbxhzy at hotmail.com