#1.首先下载结构域相互作用数据
#2.用hmmscan软件扫描所有蛋白序列的结构域
/data2/lizhaohong/bin/hmmscan --tblout 1wuzefeng.ara.hmm.out -E 1e-2 --cpu 40 /data1/wuzefeng/blast3p/Arabidopsis_thaliana.TAIR10.31.pep.all.fa_longest_isoform.fa
#3.将同一个domain的基因用python脚本进行分析
python 1hmm.parse.py >2ara_gene_inteact.txt 产生基因互作信息
#4.对结果文件进行unique
awk '$1!=$2' 2ara_gene_inteact.txt >2.5ara_gene_inteact.txt  #去除自己和自己的结构域互作关系
sort 2.5ara_gene_inteact.txt | uniq > 3ara_int.txt   # sort/uniqe

#5.用R读取后
对两列基因排序,方便训练集查找
as.data.frame(t(apply(3ara_int.order.txt,1,sort)) > 4ara_int.txt :贝叶斯整合时直接导入4ara_int.txt

