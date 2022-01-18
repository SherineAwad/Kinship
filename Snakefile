configfile: "config.yaml"

with open(config['SAMPLES']) as fp:
    SAMPLES= fp.read().splitlines()
print(SAMPLES)


rule all:
      input:
           "out.relatedness",
           "out.relatedness2",
           expand("{sample}.roh", sample=SAMPLES)
rule relatedness: 
     input: 
        vcf = expand("{sample}.vcf", sample=SAMPLES)
     log: "logs/relatedness.log" 
     benchmark: "logs/relatedness.benchmark" 
     output: 
        "out.log",
        "out.relatedness" 
     shell:
         """
         vcftools --gzvcf {input} --relatedness
         """



rule relatedness2:
     input:
        vcf = expand("{sample}.vcf", sample=SAMPLES)
     log: "logs/relatedness.log"
     benchmark: "logs/relatedness.benchmark"
     output:
        "out.log",
        "out.relatedness2"
     shell:
         """
         vcftools --gzvcf {input} --relatedness2
         """

rule ROH: 
     input:
        vcf = expand("{sample}.vcf", sample=SAMPLES)
     log: "logs/ROH.log"
     benchmark: "logs/ROH.benchmark"
     params: 
         G = config['G'],
         AF = config['AF']  
     output: expand("{sample}.roh", sample=SAMPLES) 
     shell:
        """     
	bcftools roh -G{params.G} --AF-dflt {params.AF} {input} > {output}
        """ 
