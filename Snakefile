configfile: "config.yaml"

SAMPLES = config['SAMPLES']

rule all:
      input:
           expand("{sample}.bam", sample =config['SAMPLES']),
           expand("{cohort}.vcf", cohort=config['VCF']),
           expand("{vcf}.vcf.gz", vcf=config['VCF']),
           expand("{vcf}.vcf.gz.tbi", vcf=config['VCF']),
           expand("{vcf}.relatedness",vcf=config['VCF']),
           expand("{vcf}.relatedness2", vcf=config['VCF']),
           expand("{vcf}.bed", vcf = config['VCF']),
           expand("{vcf}.kin", vcf = config['VCF']),
           expand("{sample}.kinship.txt", sample =config['SAMPLES'])

if config['PAIRED']:
    rule trim:
       input:
           r1 = "{sample}.r_1.fq.gz",
           r2 = "{sample}.r_2.fq.gz"
       output:
           "galore/{sample}.r_1_val_1.fq.gz",
           "galore/{sample}.r_2_val_2.fq.gz"
       conda: 'env/env-trim.yaml'
       shell:
           """
           mkdir -p galore
           mkdir -p fastqc
           trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore --paired {input.r1} {input.r2}
           """
    rule tosam:
        input:
            genome =expand("{genome}.fasta", genome = config['GENOME']),
            r1 = "galore/{sample}.r_1_val_1.fq.gz",
            r2 = "galore/{sample}.r_2_val_2.fq.gz"
        output:
            '{sample}.sam'
        conda: 'env/env-align.yaml'
        shell:
            "bwa mem {input.genome} {input.r1} {input.r2} > {output}"
else:
     rule trim:
       input:
           "{sample}.trimmed.fq.gz",
       output:
           "galore/{sample}_trimmed.fq.gz",
       conda: 'env/env-trim.yaml'
       shell:
           """
           mkdir -p galore
           mkdir -p fastqc
           trim_galore --gzip --retain_unpaired --trim1 --fastqc --fastqc_args "--outdir fastqc" -o galore {input}
           """
     rule tosam:
        input:
           genome =expand("{genome}.fasta", genome = config['GENOME']),
           reads = "galore/{sample}_trimmed.fq.gz"
        
        output:
            '{sample}.sam'
        conda: 'env/env-align.yaml'
        shell:
           "bwa mem {input.genome} {input.reads} > {output}"

rule sam_bam:
    input:
        "{sample}.sam"
    output:
        "{sample}.bam"
    shell:
         """
         samtools view -S -b {input} > {output}
         samtools index {input}
         """

rule tobcf: 
    input: 
       "{sample}.bam"
    output: 
       "{sample}.bcf" 
    params: 
       expand("{genome}.fasta", genome = config['GENOME'])
    shell: 
       """ 
       bcftools mpileup --fasta-ref {params} {input} -d 10000 | bcftools call -vcO v -o {output} 
       """

rule vcf:
    input:
        expand("{genome}.fasta", genome = config['GENOME']),
    params:
         I =  lambda w: " -Ou " +" ".join(expand("{sample}.bam", sample =config['SAMPLES']))
    output:
        expand("{cohort}.vcf",  cohort=config['VCF'])
    shell:
        """
        bcftools mpileup --fasta-ref {input} {params.I} -d 10000 --threads 10 | bcftools call -vcO v -o {output}
        """
 
rule bgzip:       
     input:
        "{sample}.vcf"
     output:
       "{sample}.vcf.gz"
     shell:
         """
         bgzip -c {input} > {output}
         """

rule tabix:
     input:
        "{sample}.vcf.gz"
     output:
        "{sample}.vcf.gz.tbi"
     shell:
         """
         tabix -p vcf {input}
         """

rule relatedness: 
     input: 
        "{sample}.vcf.gz"
     params: 
        "{sample}"
     output: 
        "{sample}.log",
        "{sample}.relatedness" 
     shell:
         """
         vcftools --gzvcf {input} --relatedness --out {params} 
         """


rule relatedness2:
     input:
        "{sample}.vcf.gz"
     params:
        "{sample}"
     output:
        "{sample}.log",
        "{sample}.relatedness2"
     shell:
         """
         vcftools --gzvcf {input} --relatedness2 --out {params}
         """

rule makebed: 
   input:
        "{sample}.vcf.gz"
   output:
       "{sample}.bed",
       "{sample}.fam",
       "{sample}.bim"
   params: 
      vcf = SAMPLES 
   shell:
      """
      plink --vcf {input} --make-bed --out {params} 
      """

rule kinship:
   input:
        "{sample}.bed", 
        "{sample}.fam",
        "{sample}.bim"
   params: 
       vcf = SAMPLES
   output:
       "{sample}.kin"
   shell:
      """
      king -b {input[0]} --fam {input[1]} --bim {input[2]} --related 
      king -b {input[0]} --fam {input[1]} --bim {input[2]} --kinship 
      king -b {input[0]} --fam {input[1]} --bim {input[2]} --ibdseg
      king -b {input[0]} --fam {input[1]} --bim {input[2]} --ibs 
      king -b {input[0]} --fam {input[1]} --bim {input[2]} --homog  
      mv king.kin {output}
      """
rule get_phased: 
   output: 
       "1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf", 
       "1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.idx"
   shell: 
      """
      wget wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf
      wget wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf.idx

      """  
       
rule akt_kinship: 
   input: 
       "{sample}.bcf", 
       "1000G.phase3.integrated.sites_only.no_MATCHED_REV.hg38.vcf" 
   output: 
       "{sample}.kinship.txt" 
   shell: 
       """
       akt kin -R {input[1]} -M 1 {input[0]} > {output} 
       """ 
