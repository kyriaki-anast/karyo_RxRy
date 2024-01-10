# karyo_RxRy
Script for karyotypic sex inference and detection of autosomal and sex chromosomal aneuploidies using ancient DNA.  
The script is an update to Skoglund et al. 2013 (https://github.com/pontussk/ry_compute), identifying individuals with possible 
sex karyotypes: **46XX, 46XY, 47XXX, 47XXY, 47XYY, 45X0**, as well as additional copies of chromosome 21.  
Recommended input -> **BAM files**. 

## Usage
### Whole genome sequencing data 
For shotgun sequencing data, we recommend using the default version of the script on filtered BAM files with MAPQ > 30. 
#### Example using shotgun sequencing data
``` samtools view -q 30 bamfile.bam | python2 karyo_RxRy.py ``` 
### Target enrichment data 
For genome-wide data sequenced after target enrichment with the "1240k" array, we recommend using the ```--capture``` option.
#### Example using target enrichment data
``` samtools view -q 30 bamfile.bam | python2 karyo_RxRy.py --capture ``` 
### Chromosome 21 
Using the option ```--chr21```, the output will include 2 additional columns for the R21 estimate and the Standard Error.  
Estimates higher than **0.02** indicate the presence of an additional copy of chromosome 21 and therefore a karyotype of 47X*,+21. 
#### Example command to output R21 
``` samtools view -q 30 bamfile.bam | python2 karyo_RxRy.py --chr21 ```

## Example output

| Na | Rx  | RxSE  | Ry | RySE  | Assignment  |
| :------ | :------ | :------ | :------ | :------ | :------ |
| 17243391 | 0.0297 | 4.088e-05 | 0.0024 | 1.1915e-05 | XY

  - **Na:**         Number of sequences mapping to autosomes (excluding chromosomes 13, 18, and 21)
  - **Rx:**         Number of sequences mapping to chromosome X over Na
  - **RxSE:**       Standard Error for Rx 
  - **Ry:**         Number of sequences mapping to chromosome Y over Na
  - **RySE:**       Standard Error for Ry
  - **Assignment:** Identified karyotype with possible values:
    - XX
    - XY
    - XXX
    - XXY
    - XYY
    - X0
    - Contamination
    - Consistent_with_XX
    - Consistent_with_XY

Output header can be removed using ```--noheader```. 

## Citation
Anastasiadou, K. et al. (2023). Detection of chromosomal aneuploidy in ancient genomes. Communications Biology. 10.1038/s42003-023-05642-z.

## Contacts
kyriaki.anastasiadou@crick.ac.uk  
pontus.skoglund@gmail.com
