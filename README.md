# 提取基因组中的内含子、外显子以及基因间区
## 前言
 ![image](https://user-images.githubusercontent.com/71922803/187844931-0e1ee826-be81-4f85-aa33-ad1c82502c4c.png)

DaveTang的这篇博客更新于2014年，那时转录组测序很火热。RNA-Seq往往是把reads比对回基因组或转录组，然后就是对reads进行注释，看看它们落在什么位置，其实定量也是其中一步。这里，作者就带我们看看他是如何获得外显子（exonic）、内含子（intronic）、基因间区（intergenic regions）的坐标的。
## 首先下载参考基因组注释
```
#hg38版本
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gtf.gz

zcat gencode.v33.annotation.gtf.gz | head -n6 
##description: evidence-based annotation of the human genome (GRCh38), version 33 (Ensembl 99)
##provider: GENCODE
##contact: gencode-help@ebi.ac.uk
##format: gtf
##date: 2019-12-13
chr1	HAVANA	gene	11869	14409	.	+	.	gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; hgnc_id "HGNC:37102"; havana_gene "OTTHUMG00000000961.2";
```
其中的第三列说明了feature的类型，其中会有exon、CDS、UTR等等
```
zcat gencode.v33.annotation.gtf.gz| grep -v "^#" | cut -f3 | sort | uniq -c | sort -k1rn
1377112 exon
 763824 CDS
 311868 UTR
 227912 transcript
  87848 start_codon
  80187 stop_codon
  60662 gene
    119 Selenocysteine
    ```
