# ARG analysis. The pipeline is run by ARG.smk

## I. prepare vcf files and annotations

  add_TEs_to_vcf.py
  extract_euchromatic_snps.sh
  snpEff.sub

## II. run Singer

  #Singer.sub
  #compute_trace.sub #check converge
Update 1: use the output *with_TEs.ann.haploid.vcf from haploidize_vcf.py as input to run Singer
Update 2: run Singer version 0.1.9, located in path_to_singer=/dfs7/grylee/yuhenh3/mutation_accumulation/ARG/singer-0.1.9-beta-linux-x86_64
Here is the instruction:
I want to run the entire chromosome arm (e.g. 2L). Run singer_master for continuous segments (such as 5Mb) and then use the following tool to merge them together.
For singer_master still use parameters
SINGER_PARAMS = dict(
    Ne        = "1e6",
    mu        = "5.0e-9",
    ratio     = "2.0",
)

````
path_to_singer/singer_master -Ne 1e6 -m 5.0e-9 -ratio 2.0 -vcf prefix_of_vcf_file -output prefix_of_output_file_0 -start 0 -end 5e6
path_to_singer/singer_master -Ne 1e6 -m 5.0e-9 -ratio 2.0 -vcf prefix_of_vcf_file -output prefix_of_output_file_1 -start 5e6 -end 10e6
path_to_singer/singer_master -Ne 1e6 -m 5.0e-9 -ratio 2.0 -vcf prefix_of_vcf_file -output prefix_of_output_file_2 -start 10e6 -end 15e6
......
````
`python merge_ARG.py --file_table sub_file_table_file --output merged_ARG_filename`
The sub_file_table_file specifies how the inferred ARG should be pieced together with an example like this:
````
prefix_of_output_file_0_nodes_i.txt prefix_of_output_file_0_branches_i.txt prefix_of_output_file_0_muts_i.txt 0 
prefix_of_output_file_1_nodes_i.txt prefix_of_output_file_1_branches_i.txt prefix_of_output_file_1_muts_i.txt 5000000
prefix_of_output_file_2_nodes_i.txt prefix_of_output_file_2_branches_i.txt prefix_of_output_file_2_muts_i.txt 10000000
......
````

where the middle index is the index of the genomic window on which parallelization is performed, and i is the index of the ARG replicate (MCMC sample). Note that the genomic windows do not have to be all adjacent. For example, when there is a centromeric region in the middle of the chromosome, it can be skipped by only selecting windows in the left and right arms. This will be reflected in the last column of the sub_file_table_file, which indicates the starting position of each genomic window.

By doing the merging operation, you will get the ARG samples with the length of the entire chromosome rather than for each individual genomic window, which might be convenient in certain scenarios.



## III. redate
Update 3: the new outputs from Singer will run remove_TE_from_trees.py, only remove the TEs, keep all SNP mutations. Then redate them with Polegon

  remove_TE_from_trees.py #remove TEs from the tree
  Polegon.sub


## IV. analysis
Update 4: only use rule mutation_rate_weight_on_TE_trees, for TE, strong SNP mutation and neutral SNP mutation.


  Analysis.sub
