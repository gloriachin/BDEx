# Usage:

$ from genesets import get_geneset

## get genesets from msidgdb
Example
$ BDEx.get_geneset("msigdb","HALLMARK_APOPTOSIS")

## get datasets from Depmap/CCLE etc
Example:
$ BDEx.get_data("Depmap_21Q4","sample_info")

The following lists are currently supported.
Depmap_21Q4,CCLE_gene_cn<br />
Depmap_21Q4,CCLE_expression<br />
Depmap_21Q4,CCLE_mutation<br />
Depmap_21Q4,sample_info<br />
Depmap_21Q4,CRISPR_gene_effect<br />
Depmap_21Q4,CRISPR_gene_dependency<br />
Depmap_21Q4,CRISPR_common_essentials<br />
Depmap_21Q4,common_essentials<br />
DEMETER2_V6,sample_info<br />
DEMETER2_V6,D2_DRIVE_gene_dep_scores<br />
DEMETER2_V6,D2_combined_gene_dep_scores<br />
