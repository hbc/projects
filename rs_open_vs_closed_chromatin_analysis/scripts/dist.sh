samples="adenoma1
adenoma2
adenoma3
adenoma4
adenoma5
carcinoma1
carcinoma2
carcinoma3
carcinoma4
carcinoma5"



for sample in $samples
do
IntervalStats -D -q ~/bin/intervalstats/test_dists.bed  -r open_${sample}_top1000.bed -d ./promoters/promotes.500.5000.bed -o ${sample}.promoters.dist
done

for sample in $samples
do
IntervalStats -D -q ~/bin/intervalstats/test_dists.bed  -r open_${sample}_top1000.bed -d ./wg/gold.hg18.genome.bed -o ${sample}.wg.dist
done