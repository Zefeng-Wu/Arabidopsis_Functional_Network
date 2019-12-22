#step1:
下载5个基因组蛋白数据库
#step2:
三个物种全蛋白序列比对到各基因组数据,得到blast_out结果
#step3:
用5个python脚本分别处理blast结果,得到phylogenetic_profile目录,存放gene-species-相似度
#step4:
用R脚本聚类去冗余物种,得到部分数据,放到gene2species目录



