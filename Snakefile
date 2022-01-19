configfile: "config.yaml"

VCF = config['VCF']

rule all:
      input:
           expand("{sample}.vcf.gz", sample=VCF),
           expand("{sample}.vcf.gz.tbi", sample=VCF),
           expand("{sample}.relatedness",sample=VCF),
           expand("{sample}.relatedness2", sample=VCF),
           expand("{sample}.bed", sample = VCF),
           expand("{sample}.kin", sample = VCF)
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
      vcf = VCF 
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
       vcf = VCF
   output:
       "{sample}.kin"
   shell:
      """
      king -b {input[0]} --fam {input[1]} --bim {input[2]} --related --kinship, --ibdseg, --ibs, --homog  
      mv king.kin {output}
      """
