# BSA-cotton
  NGS Bulk Segregant Analysis in Gossypium hirsutum L
 #先对基因组文件建立索引
bwa index AD1_TM-1_UTX.2.0.fa
#样本中双亲为父本：Maxxa，母本：TX2094-SZ，子代极端池分别为Flower1，Flower3
#对测序数据质控-fastp原始数据进行一波过滤
mkdir -p 01-clean-data
fastp -i data/Flower1/Flower1_1.fq.gz -I data/Flower1/Flower1_2.fq.gz -o 01-clean-data/Flower1_1.fq.gz -O 01-clean-data/Flower1_2.fq.gz
fastp -i data/Flower3/Flower3_1.fq.gz -I data/Flower3/Flower3_2.fq.gz -o 01-clean-data/Flower3_1.fq.gz -O 01-clean-data/Flower3_2.fq.gz
fastp -i data/maxxa/Maxxa_1.fq.gz -I data/maxxa/Maxxa_2.fq.gz -o 01-clean-data/Maxxa_1.fq.gz -O 01-clean-data/Maxxa_2.fq.gz
fastp -i data/TX2094-SZ/TX2094-SZ_raw_R1.fq.gz -I data/TX2094-SZ/TX2094-SZ_raw_R2.fq.gz -o 01-clean-data/TX2094-SZ_raw_R1.fq.gz -O 01-clean-data/TX2094-SZ_raw_R2.fq.gz
#数据质控，输出HTML文件
fastqc Flower2/Flower2_1.fq.gz Flower2/Flower2_2.fq.gz -o fQreport/ -t 4
#对过滤后的文件进行序列比对，回帖到参考基因组使用BWA
#给变量赋值
NUMBER_BWA_THREADS=8
NUMBER_SAMTOOLS_THREADS=4
REFERENCE="/public/agis/huguanjing_group/hushiyin/BSA/ref/AD1_TM-1_UTX.3.0.fa"
CLEAN_DATA_DIR="01-clean-data"
ALIGN_DIR="02-read-align"
mkdir -p 02-read-align
sample="Flower1 Flower3 Maxxa TX2094-SZ"
R1="${CLEAN_DATA_DIR}/${sample}_1.fq.gz"
R2="${CLEAN_DATA_DIR}/${sample}_2.fq.gz"
echo "Processing sample: ${sample}"
# BWA比对 → samtools排序 → 索引
bwa mem -M -R "@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA" \
        -t ${NUMBER_BWA_THREADS} ${REFERENCE} ${R1} ${R2} \
        | samtools sort -@ ${NUMBER_SAMTOOLS_THREADS} \
            -o "${ALIGN_DIR}/${sample}_sort.bam" -
    #建立BAM索引
    samtools index -@ ${NUMBER_SAMTOOLS_THREADS} \
        "${ALIGN_DIR}/${sample}_sort.bam"
#这一步可以进行质控，筛选Q>30的单一比对位点获得unique.bam,后续进行去重复与变异检测#
amtools view -@ 4 -b -h -q 30 -F 0x904 \
    Flower1_sort.bam | \
    samtools sort -@ 4 -o Flower1_unique_sort.bam

samtools index -@ 4 Flower1_unique_sort.bam 
#去除相同位置且序列一模一样的reads，去除重复序列，拥有相同位置信息，且序列也一模一样的read，通常是由PCR扩增引起
sambamba markdup -r -t 8 02-read-align/Flower1_sort.bam 02-read-align/Flower1_markdup.bam
sambamba markdup -r -t 8 02-read-align/Flower3_sort.bam 02-read-align/Flower3_markdup.bam
sambamba markdup -r -t 8 02-read-align/Maxxa_sort.bam 02-read-align/Maxxa_markdup.bam
sambamba markdup -r -t 8 02-read-align/TX2094-SZ_sort.bam 02-read-align/TX2094-SZ_markdup.bam
# 使用GATK进行SNP calling
for sample in Flower1 Flower3 Maxxa TX2094-SZ; do
  gatk HaplotypeCaller \
    -R ../ref/AD1_TM-1_UTX.3.0.fa \
    -I ${sample}_markdup.bam \
    -O ../03-variants/gvcf/${sample}.g.vcf.gz \
    --emit-ref-confidence GVCF \               #生成 GVCF 格式（用于 joint genotyping）
    --dont-use-soft-clipped-bases \            #忽略 soft-clipped reads，减少假阳性
    -stand-call-conf 30 \                      #设定变异调用置信度阈值
    --native-pair-hmm-threads 8 \              #启用 8 线程以加速计算
    > ../logs/${sample}.log 2>&1
done
#Identify variation进行SNP的识别变异
#首先对不同类型的混池的四个VCF文件进行合并
gatk CombineGVCFs \
  -R ../ref/TM1.UTX.3.0.fa \
  -V gvcf/Maxxa.g.vcf.gz \
  -V gvcf/TX2094-SZ.g.vcf.gz \
  -V gvcf/Flower5.g.vcf.gz \
  -V gvcf/Flower7.g.vcf.gz \
  -O gatk13MT.g.vcf.gz
#Identify 与select SNP
gatk GenotypeGVCFs \
    -R ../ref/TM1.UTX.3.0.fa \
    -V gatk13MT.g.vcf.gz \
    -O gatk13MT.vcf.gz
    #只保留 SNP
gatk SelectVariants \
    -R ../ref/TM1.UTX.3.0.fa \
    -V gatk13MT.vcf.gz \
    --select-type-to-include SNP \
    -O gatk13MT.snp.vcf.gz
#Filter过滤SNP
gatk VariantFiltration \
    -R ../ref/TM1.UTX.3.0.fa \
    -V gatk13MT.snp.vcf.gz \
    --missing-values-evaluate-as-failing \      #缺失值时视为不合格
    --cluster-window-size 10 \    #控制聚类过滤（防止过度聚集的突变）在10 bp 范围内检测
    --cluster-size 3 \           #给定窗口中发现 ≥3 个变异，则认为这些变异是可疑的
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \      #(Quality by Depth）(Mapping Quality)(Fisher Strand Bias)(Strand Odds Ratio)变异等位基因和参考等位基因之间的比对质量差异，变异位点在 reads 中的位置差异
    --filter-name "Filter" \
    -O gatk13MT.snp.filter.vcf.gz
