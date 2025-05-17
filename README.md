# Mitogenomes


# 1) mapeando los raw reads con el metagenoma

#!/bin/bash
#SBATCH --job-name=mitogenome_boitata
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --tasks-per-node=15 # Number of cores per node
#SBATCH --time=40:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgarcia@amnh.org
#SBATCH --output=slurm-%j-%x.out
#SBATCH --error=slurm-%j-%x.err


source ~/.bash_profile
module load BWA/bwa-0.7.17

#escritorio de trabajo
cd /home/dgarcia/mendel-nas1/short_reads/mitogenomes

bwa mem -t 14 /home/dgarcia/mendel-nas1/short_reads/mitogenomes/Erythrolapmus.fasta \
Helicops_boitata_CHUNB_83210_R1_clean.fastq.gz \
Helicops_boitata_CHUNB_83210_R2_clean.fastq.gz > boitata_mapped.sam
 
 
 hagmanni_mapped.sam
 
 # 2) transformando el .sam file a un .bam file 
 
#!/bin/bash
#SBATCH --job-name=sam_boitata
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --tasks-per-node=15 # Number of cores per node
#SBATCH --time=50:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgarcia@amnh.org
#SBATCH --output=slurm-%j-%x.out
#SBATCH --error=slurm-%j-%x.err

source ~/.bash_profile
module load Samtools/samtools-1.11

# Cambiar al directorio de trabajo
cd /home/dgarcia/mendel-nas1/short_reads/mitogenomes


# 1. Convertir SAM a BAM
samtools view -@ 14 -bS boitata_mapped.sam > boitata_mapped.bam

# 2. Ordenar el BAM
samtools sort -@ 14 -o boitata_mapped.sorted.bam boitata_mapped.bam

# 3. Filtrar por MAPQ ≥ 30
samtools view -@ 14 -b -q 30 boitata_mapped.sorted.bam > boitata_mapped.filtered.bam

# 4. Indexar el BAM filtrado (opcional pero recomendado)
samtools index boitata_mapped.filtered.bam
 
 
 
 # 3) estamos haciendo un genoma consenso:
 
#!/bin/bash
#SBATCH --job-name=sam_boitata
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --tasks-per-node=15 # Number of cores per node
#SBATCH --time=50:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dgarcia@amnh.org
#SBATCH --output=slurm-%j-%x.out
#SBATCH --error=slurm-%j-%x.err
 
 #escritorio de trabajo
cd /home/dgarcia/mendel-nas1/short_reads/mitogenomes

source ~/.bash_profile
conda activate bcftools

REF="/home/dgarcia/mendel-nas1/short_reads/mitogenomes/Erythrolapmus.fasta"
BAM="boitata_mapped.filtered.bam"
VCF="boitata.variants.vcf.gz"
CONSENSO="boitata.mapped.consensus.fasta"

bcftools mpileup -f "$REF" "$BAM" -Ou \
  | bcftools call -mv -Oz -o "$VCF"

bcftools index "$VCF"

bcftools consensus -f "$REF" "$VCF" > "$CONSENSO"

echo " Secuencia consenso creada: $CONSENSO"

# Otros comando para mirar la cualidad del genoma:

metricas para pillar que tan bien quedaron alineados los mitogenomas 

a) Samtools flagstat mapped.filtered.bam

b) Cobertura por posición:

samtools depth -a mapped.filtered.bam > cobertura.txt
samtools depth -a Helicops_hagmanni_CHUNB_83215_mapped.filtered.bam > hagmanni.cobertura.txt

# Luego calculamos la cobertura promedio

awk '{sum+=$3} END { print "Cobertura promedio:",sum/NR }' cobertura.txt
awk '{sum+=$3} END { print "Cobertura promedio:",sum/NR }' cobertura.txt


c) samtools stats mapped.filtered.bam > mapeo.stats

d) plotear la cobertura (Rstudio)
cov <- read.table("cobertura.txt")
plot(cov$V2, cov$V3, type="l", xlab="Position", ylab="Coverage", main="Coverage per position")


