# ARG analysis. The pipline is run by ARG.smk

## I. prepare vcf files and annotations

  add_TEs_to_vcf.py
  extract_euchromatic_snps.sh
  snpEff.sub

## II. run Singer

  Singer.sub
  #compute_trace.sub #check converge

## III. redate

  remove_TE_from_trees.py #remove TEs from the tree
  Polegon.sub

## IV. analysis

  Analysis.sub
