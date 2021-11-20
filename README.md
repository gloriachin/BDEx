# Usage:

$ from genesets import get_geneset

## get genesets from msidgdb
Example
$ BDEx.get_geneset("msigdb","HALLMARK_APOPTOSIS")

## get datasets from Depmap/CCLE etc
Example:
$ BDEx.get_data("Depmap_21Q4","sample_info")

The following lists are currently supported.
Depmap_21Q4,CCLE_gene_cn
Depmap_21Q4,CCLE_expression
Depmap_21Q4,CCLE_mutation
Depmap_21Q4,sample_info
Depmap_21Q4,CRISPR_gene_effect
Depmap_21Q4,CRISPR_gene_dependency
Depmap_21Q4,CRISPR_common_essentials
Depmap_21Q4,common_essentials
DEMETER2_V6,sample_info
DEMETER2_V6,D2_DRIVE_gene_dep_scores
DEMETER2_V6,D2_combined_gene_dep_scores
