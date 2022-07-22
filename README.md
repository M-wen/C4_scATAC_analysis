# C4_scATAC_analysis
# 数据预处理

## DNBelabC4 scATAC数据结构

![](image/reads_structure.png)

### 比对

使用chromap进行比对 chromap[**参考文献**](https://www.nature.com/articles/s41467-021-26865-w)

chromap 的优点

-   比对速度较bwa mem2 提升了10倍

-   比对准确性与bwa mem2相当

![chromap 比对原理](image)

### beads合并

使用d2c根据beads中fragment中的相似性进行beads 合并，[**参考文献**](https://www.nature.com/articles/s41467-020-14667-5)

-   液滴中出现捕捉多个beads一个细胞的情况

![一个液滴 多beads或多个barcode情况](image/drop_%E5%A4%9Abeads.png){width="594"}

-   通过计算两个beads之间相似性，推测beads来源于同一个液滴

![计算两个barcode之间相似性](image/bap%E5%8E%9F%E7%90%86.png){width="698"}

### peak calling

![](image/peak_calling_MACS2.png)
