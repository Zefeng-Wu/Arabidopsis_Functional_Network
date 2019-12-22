#扫描植物基因启动子转录因子结合位点TFBS

#从jaspar网站下载植物转录因子motif,core data中的fpm文件:http://jaspar.genereg.net/html/DOWNLOAD/JASPAR_CORE/pfm/nonredundant/pfm_plants.txt/ 放入目录TF_binding_motif-database
sed -e 's/^[ACGT]//g' -e 's/\[//g' -e 's/\]//g' pfm_plants.pfm >1plants.pfam #转换格式
jaspar2meme -bundle  1plants.pfm >plants.meme #转成meme格式
fimo plants.meme ../../TAIR10_data/2ara_promoters.fa #扫描启动子区

