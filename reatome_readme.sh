reactome下载拟南芥蛋白互作数据
http://www.reactome.org/download/current/
awk 'BEGIN{FS="\t"}{print $2,$5}' arabidopsis_thaliana.interactions.txt | sort | uniq | awk '$1!=$2' > ara.reactome

1.提取第2列和第5列数据
2.把一对多用脚本分开
