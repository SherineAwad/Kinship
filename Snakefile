configfile: "config.yaml"

with open(config['SAMPLES']) as fp:
    SAMPLES= fp.read().splitlines()
print(SAMPLES)


rule all:
      input:
           expand("{sample}.relatedness",sample=SAMPLES),
           expand("{sample}.relatedness2", sample=SAMPLES),
           expand("{sample}.roh", sample=SAMPLES), 
           expand("{sample}.vcf.gz", sample=SAMPLES), 
           expand("{sample}.vcf.gz.tbi", sample=SAMPLES), 
           expand("{all}.vcf.gz", all = config['ALL']), 
           expand("{sample}.stats", sample=SAMPLES)
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
        "{sample}.vcf"
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
        "{sample}.vcf"
     params:
        "{sample}"
     output:
        "{sample}.log",
        "{sample}.relatedness2"
     shell:
         """
         vcftools --gzvcf {input} --relatedness2 --out {params}
         """


rule ROH: 
     input:
        "{sample}.vcf"
     params: 
         G = config['G'],
         AF = config['AF']  
     output: "{sample}.roh" 
     shell:
        """     
	bcftools roh -G{params.G} --AF-dflt {params.AF} {input} > {output}
        """
rule merge: 
    input: 
        vcf = expand("{sample}.vcf.gz", sample=SAMPLES) 
    output:
        expand("{all}.vcf.gz", all = config['ALL'])
    shell: 
        """ 
        vcf-merge  {input} | bgzip -c > {output} 
        """ 

rule stats: 
   input:
        "{sample}.vcf.gz"
   output:
        "{sample}.stats"
   shell: 
      """ 
      vcf-stats {input} > {output} 
      """  
