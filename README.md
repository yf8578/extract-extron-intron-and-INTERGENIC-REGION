# 提取基因组中的内含子、外显子以及基因间区
文章来自：https://www.jieandze1314.com/post/cnposts/165/ ：仅供本人学习参考
## 前言

<div align=center>
<img src="https://user-images.githubusercontent.com/71922803/187844931-0e1ee826-be81-4f85-aa33-ad1c82502c4c.png" width="450">
</div>


DaveTang的这篇博客更新于2014年，那时转录组测序很火热。RNA-Seq往往是把reads比对回基因组或转录组，然后就是对reads进行注释，看看它们落在什么位置，其实定量也是其中一步。这里，作者就带我们看看他是如何获得外显子（exonic）、内含子（intronic）、基因间区（intergenic regions）的坐标的。
## 首先下载参考基因组注释
```
#hg38版本
wget -c ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.annotation..gz

zcat gencode.v41.annotation.gtf.gz | head -n6 
##description: evidence-based annotation of the human genome (GRCh38), version 41 (Ensembl 107)
##provider: GENCODE
##contact: gencode-help@ebi.ac.uk
##format: gtf
##date: 2022-05-12
chr1    HAVANA  gene    11869   14409   .       +       .       gene_id "ENSG00000223972.5"; gene_type "transcribed_unprocessed_pseudogene"; gene_name "DDX11L1"; level 2; hgnc_id "HGNC:37102"; havana_gene "OTTHUMG00000000961.2";
```
其中的第三列说明了feature的类型，其中会有exon、CDS、UTR等等
```
zcat gencode.v44.annotation.gtf.gz| grep -v "^#" | cut -f3 | sort | uniq -c | sort -k1rn
1625321 exon
 871490 CDS
 377983 UTR
 251236 transcript
  97009 start_codon
  90749 stop_codon
  61852 gene
    119 Selenocysteine
```
它们之间的关系，可以看看：http://blog.sciencenet.cn/blog-1271266-792469.html 以及 http://www.dxy.cn/bbs/topic/36728037
UTR（UntranslatedRegions)即非翻译区，是信使RNA（mRNA）分子两端的非编码片段，`UTR在DNA序列中算是外显子`;5’-UTR从mRNA起点的甲基化鸟嘌呤核苷酸帽延伸至AUG起始密码子，3’-UTR从编码区末端的终止密码子延伸至多聚A尾巴（Poly-A）的末端；
基因间区是每个基因之间的间隔序列，不属于外显子也不属于内含子，它是non coding；

<div align=center>
<img src="https://user-images.githubusercontent.com/71922803/187847850-3b5d9a82-d988-49af-b303-775c5456a26d.png" width="450">
</div>

>【图片中描述的过程：基因经过转录形成 pre-mRNA，这里面包含着内含子和外显子（5端是一个外显子，但是这段外显子不仅包含CDS，还包含5’ UTR；3端也是以外显子结束，但它也是不仅包含CDS，还包含3’ UTR），经过剪接形成成熟mRNA,内含子已减掉。如果抛开后来加上去的cap和poly A的话，这时全是外显子，但是不全是CDS，因为只有中间的那部分以起始密码子开始、以终止密码子结束的片段才是CDS，只有这部分才会被翻译成蛋白质】

> **转录本（transcript)**是基因通过转录形成的一种或多种可供编码蛋白质的成熟的mRNA。一个基因由于可变剪切现象，可能产生多个转录本。【基因转录的过程，首先是形成前体mRNA，通过剪切内含子连接外显子，5’端加帽及3’端加尾之后形成成熟的mRNA。可变剪切就是：在剪切的过程中可能会剪切掉外显子，也有可能保留部分内含子，因此会形成多种mRNA即多个转录本】，而我们平时研究某个基因的功能，实际是研究它的一个转录本编码的蛋白的功能。一般情况下它的不同的转录本分布在不同类型的细胞中，当然也有可能多种转录本同时存在于某一细胞中。
>编码区（Coding region，CDS） 是mRNA序列中编码蛋白质的那部分序列
## 获得外显子、内含子以及基因间区
### 外显子的获得 => mergeBed
从GTF文件中其实已经能找到exon的踪迹，但是同一个基因下的各个exon可能会有交集，于是这里作者在统计时，将存在重叠的外显子融合在一起（这里只是采用最简单的方法进行探索，但不是最准确的）
<div align=center>
<img src="https://user-images.githubusercontent.com/71922803/187849154-f6bde10d-d0da-403e-a1af-4d6d5dcdf0c8.png" width="450">
</div>

使用`bedtools`的`mergeBed`功能，但前提是先使用`sortBed`将坐标进行排序

```
# 首先挑出exon的起始和终止位点
zcat gencode.v33.annotation.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5}' | head -n5
chr1	11868	12227
chr1	12612	12721
chr1	13220	14409
chr1	12009	12057
chr1	12178	12227
# 然后进行排序，接着进行merge
zcat gencode.v33.annotation.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="exon" {print $1,$4-1,$5}' | sortBed | mergeBed -i - | gzip >gencode.v33.exon.annotation.gtf.gz

zcat gencode.v33.exon.annotation.bed.gz | head -n5
chr1	11868	12227
chr1	12612	12721
chr1	12974	13052
chr1	13220	14501
chr1	15004	15038
```
>小tip： 上面👆的awk命令中出现了一个{print $1,$4-1,$5}，其中$4是指的GTF的第4列，即起始位点。它要减一是因为要从GTF变为BED，而它们的坐标系统不一致：GTF has 1-based start coordinates, and BED has 0-based start coordinates

### 内含子的获得 => subtractBed
<div align=center>
<img src="https://user-images.githubusercontent.com/71922803/187851318-9ffa4ff4-8fd2-4657-b0d1-84814248f42b.png" width="450">
</div>

```
# 选择GTF中gene的feature，然后排除merged exon
zcat gencode.v33.annotation.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' | sortBed | subtractBed -a stdin -b gencode.v33.exon.annotation.bed.gz | gzip >gencode.v33.intron.annotation.bed.gz

zcat gencode.v33.intron.annotation.bed.gz | head -n5
# 结果可以和上面👆的exon结果对照
chr1	12227	12612 
chr1	12721	12974
chr1	13052	13220
chr1	14501	15004
chr1	15038	15795
```
### 基因间区的获得 => complementBed

<div align=center>
<img src="https://user-images.githubusercontent.com/71922803/187851565-67827656-4875-4f35-bdc7-bb04ecdd724c.png" width="450">
</div>

需要先获得整个基因组各个染色体的长度，然后从中排除掉gene

```
# 各个染色体的长度(在bedtools complementBed的帮助文档中)
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e         "select chrom, size from hg38.chromInfo" >hg38.genome
# 去掉hg38.genome第一行注释：chrom	size
 sed -i '1d' hg38.genome
# 排除gene
zcat gencode.v33.annotation.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' | sortBed | complementBed -i stdin -g hg38.genome | gzip >gencode.v33.intergenic.annotation.bed.gz
```
**但这里报了个错：Error: Sorted input specified, but the file stdin has the following record with a different sort order than the genomeFile hg38.genome**
**实际上问题就是因为我们的sort bed后的chrxx不是按照数字排序的（http://asearchforsolutions.blogspot.com/2018/11/error-sorted-input-specified-but-file.html**
```
zcat gencode.v33.annotation.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' | sortBed | cut -f1 | uniq
chr1
chr10
chr11
chr12
chr13
chr14
chr15
chr16
chr17
chr18
# 而这里的hg38.genome却是拍好顺序的
cat hg38.genome | cut -f1 | uniq | head
chr1
chr2
chr3
chr4
chr5
chr6
chr7
chrX
chr8
chr9
# 因此我们“恃强凌弱”，将hg38.genome改成和sortBed后一致的（相比之下hg38.genome更容易修改），也不按照字母排序
cat hg38.genome| sort -k1.4 >hg38-2.genome
head hg38-2.genome
chr1
chr10
chr10_GL383545v1_alt
chr10_GL383546v1_alt
chr10_KI270824v1_alt
chr10_KI270825v1_alt
chr10_KN196480v1_fix
chr10_KN538365v1_fix
chr10_KN538366v1_fix
chr10_KN538367v1_fix
```
于是，继续进行代码，就没有问题啦：
```
zcat gencode.v33.annotation.gtf.gz | awk 'BEGIN{OFS="\t";} $3=="gene" {print $1,$4-1,$5}' | sortBed | complementBed -i stdin -g hg38-2.genome | gzip >gencode.v33.intergenic.annotation.bed.gz

zcat gencode.v33.intergenic.annotation.bed.gz | head -n5
chr1	0	11868
chr1	31109	34553
chr1	36081	52472
chr1	53312	57597
chr1	64116	65418
```
**到此便成功提取基因组中的内含子、外显子以及基因区间**
## 获得对应的.fa文件
1.下载人类基因组序列
2.使用bedtools
```
Tool:    bedtools getfasta (aka fastaFromBed)
Version: v2.30.0
Summary: Extract DNA sequences from a fasta file based on feature coordinates.

Usage:   bedtools getfasta [OPTIONS] -fi <fasta> -bed <bed/gff/vcf>

Options:
        -fi             Input FASTA file
        -fo             Output file (opt., default is STDOUT
        -bed            BED/GFF/VCF file of ranges to extract from -fi
        -name           Use the name field and coordinates for the FASTA header
        -name+          (deprecated) Use the name field and coordinates for the FASTA header
        -nameOnly       Use the name field for the FASTA header
        -split          Given BED12 fmt., extract and concatenate the sequences
                        from the BED "blocks" (e.g., exons)
        -tab            Write output in TAB delimited format.
        -bedOut         Report extract sequences in a tab-delimited BED format instead of in FASTA format.
                        - Default is FASTA format.
        -s              Force strandedness. If the feature occupies the antisense,
                        strand, the sequence will be reverse complemented.
                        - By default, strand information is ignored.
        -fullHeader     Use full fasta header.
                        - By default, only the word before the first space or tab
                        is used.
        -rna    The FASTA is RNA not DNA. Reverse complementation handled accordingly.
```
## 在提取基因间区序列过程中出现如下问题

```
WARNING. chromosome (chr10_GL383545v1_alt) was not found in the FASTA file. Skipping.
WARNING. chromosome (chr10_GL383546v1_alt) was not found in the FASTA file. Skipping.
WARNING. chromosome (chr10_KI270824v1_alt) was not found in the FASTA file. Skipping.
WARNING. chromosome (chr10_KI270825v1_alt) was not found in the FASTA file. Skipping.
WARNING. chromosome (chr10_KN196480v1_fix) was not found in the FASTA file. Skipping.
WARNING. chromosome (chr10_KN538365v1_fix) was not found in the FASTA file. Skipping.
WARNING. chromosome (chr10_KN538366v1_fix) was not found in the FASTA file. Skipping.
```

![image](https://user-images.githubusercontent.com/71922803/187893910-7c99ba00-0068-4888-a856-9f3e309707b2.png)

网上找到其他教程说是基因碎片，影响应该不大
可以根据上述图片将其剔除
