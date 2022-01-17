configfile: "config.yaml"

with open(config['SAMPLES']) as fp:
    SAMPLES= fp.read().splitlines()
print(SAMPLES)


rule all:
      input:
           "out.relatedness",
           "out.relatedness2"

rule relatedness: 
     input: 
        vcf = expand("{COHORT}.vcf", COHORT=SAMPLES)
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
        vcf = expand("{COHORT}.vcf", COHORT=SAMPLES)
     log: "logs/relatedness.log"
     benchmark: "logs/relatedness.benchmark"
     output:
        "out.log",
        "out.relatedness2"
     shell:
         """
         vcftools --gzvcf {input} --relatedness2
         """


