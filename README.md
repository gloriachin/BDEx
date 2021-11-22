# Usage:

To check the latest version of the package, please check https://test.pypi.org/project/BDEx/

Intall using pip:<br />
pip install -i https://test.pypi.org/simple/ BDEx <br />

Run the following examples in python3 environment.<br />

$ from BDEx import get_geneset_data

# Data accessibility module (get_geneset_data) <br />
## get genesets from msidgdb
Example
$ get_geneset_data.get_geneset("msigdb","HALLMARK_APOPTOSIS")

## get datasets from public datasets <br />
Example:
$ get_geneset_data.get_data("Depmap_21Q4","sample_info")

The following lists are currently supported.<br />
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

# Notes:
This datasets are identical with the orignal data resource. Users should acknolwege the data generators and use the data according to the owner's guidence. 

Link to the data resources:<br />
https://depmap.org/portal/<br />
https://www.gsea-msigdb.org/gsea/msigdb/<br />

# Analysis module (ana)

Currently provided modules: 

## 1, spearman correlation by given a dictionary type.<br />
$ ana.measure_cor_spearman_batch(); see example below.<br />

It allows a faster computation especially with large number of subjects, such as the correlation for all genes. 
## 2, check the similarity of genes or sample similarity from the co-expression results.<br />
$ ana.check_gene_groups_from_correlation(); see example below.<br />

Example:<br />
geneset = get_geneset_data.get_geneset("msigdb","HALLMARK_OXIDATIVE_PHOSPHORYLATION")<br />
df =CCLE_expr<br />

            col_names = []
            if "Unnamed: 0" in df.columns:
                for col in df.columns:
                    if col != "Unnamed: 0":
                        col_names.append(col.split(" ")[0])

                df.drop(columns=df.columns[0], 
                        axis=1, 
                        inplace=True)
            else:
                for col in df.columns:
                    if col != "Unnamed: 0":
                        col_names.append(col.split(" ")[0])

            df.columns = col_names
            geneset = list(set(geneset).intersection(set(col_names)))
            x = df.loc[:,geneset]

dic_gene_epxr = ana.pd2dic(x)<br />
result = ana.measure_cor_spearman_batch(dic_gene_epxr,geneset,geneset)<br />
[g1,g2,g3] = ana.check_gene_groups_from_correlation(result,0.3)<br />


# knowledge graph module (kg)
### Query which genes are associated with sensitivity or resistance to a drug 
kg.drug_response_KG_by_drug("JQ1","biolink:associated with sensitivity to")
kg.drug_response_KG_by_drug("Cisplatin","biolink:associated with resistance to")

### Query which drugs sensitivities are affected by the gene

drug_response_KG_by_gene("EGFR","biolink:associated with sensitivity to")
drug_response_KG_by_gene("KRAS","biolink:associated with resistance to")