import pandas as pd
import requests
import json
from scipy import stats 
import statsmodels.stats.multitest as multi
import numpy as np

def parse_result(respText):
    object_list = []
    object_name_list = []
    subject_list = []
    subject_name_list = []
    rho_list = []
    tissue_list = []
    predicates_list = []
    attribute_list =[]
    tumor_type_list = []
    method_list =[]
    modifer_list=[]
    dict_nodes = respText['message']['knowledge_graph']['nodes']
    
    for k in respText['message']['knowledge_graph']['edges'].keys():
        object_ = respText['message']['knowledge_graph']['edges'][k]['object']
        object_name = dict_nodes[object_]['name']
        
        subject_ = respText['message']['knowledge_graph']['edges'][k]['subject']
        subject_name =  dict_nodes[subject_]
        predicates = respText['message']['knowledge_graph']['edges'][k]['predicate']
        attribute = respText['message']['knowledge_graph']['edges'][k]['edge_attributes']
        
        modifier = attribute['Edge_attribute_Subject_Modifier']
        modifer_list.append(modifier)
        method = attribute['Edge_attribute_method']
        method_list.append(method)
        
        tumor_type = attribute['Edge_attribute_sample_orign']
        tumor_type_list.append(tumor_type)
        
        object_list.append(object_)
        object_name_list.append(object_name)
        
        subject_list.append(subject_)
        subject_name_list.append(subject_name)
        
        predicates_list.append(predicates)
        attribute_list.append(attribute)
        
        #rho_list.append(rho)
        #tissue_list.append(tissue_type)
    
    result = pd.DataFrame({"Subject": subject_list,"Subject_name":subject_name_list,
                           "Object":object_list, "Object_name":object_name_list , 
                           "predicate":predicates_list,"tumor_type":tumor_type_list,"Modifer":modifer_list})
    return(result)

def drug_response_KG_by_gene(Gene,  action):
    #Gene are the symbols
    #actions can be the following options: 
    # biolink:associated with sensitivity to
    # biolink:associated with resistance to

    json_query={'message': 
                    {'query_graph': 
                        {'edges': {'e00': {'subject': 'n00',
                                            'object': 'n01',
                                            'predicates': [action],

                                          }},
                         'nodes': {'n00': {'categories': ['biolink:Gene'], 
                                            'ids': ['Symbol:'+Gene]},
                                   'n01': {'categories': ['biolink:SmallMolecule']}}
                    } 
                       }}


    url = "http://35.233.133.157:5000/BigGIM/DrugResponse"
    resp = requests.request("POST", url, json=json_query)   
    respText = json.loads(resp.text)
    result = parse_result(respText)
    #result_1 = result[["Object_name","predicate","tumor_type","Modifer"]]

    return(result)

def drug_response_KG_by_drug(drug,  action):
    #Gene are the symbols
    #actions can be the following options: 
    # biolink:associated with sensitivity to
    # biolink:associated with resistance to

    json_query={'message': 
                    {'query_graph': 
                        {'edges': {'e00': {'subject': 'n00',
                                            'object': 'n01',
                                            'predicates': [action],

                                          }},
                         'nodes': {'n00': {'categories': ['biolink:Gene']},
                                   'n01': {'categories': ['biolink:SmallMolecule'],
                                           'ids': [drug]}}
                    } 
                       }}

    url = "http://35.233.133.157:5000/BigGIM/DrugResponse"
    resp = requests.request("POST", url, json=json_query)   
    respText = json.loads(resp.text)
    result = parse_result(respText)
    #result_1 = result[["Object_name","predicate","tumor_type","Modifer"]]

    return(result)

#Step 3: narrow down potential gene essentiality
def select_all_potential_essential_genes(CRISPR_gene_dependency, threshod_probability):
    # Example: select_all_potential_essential_genes(CRISPR_gene_dependency, 0.5)
    # 0.5 represents 50% of chance this gene shows essentiality in one cell line. 

    df =CRISPR_gene_dependency
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
    df.index = df['DepMap_ID']
    df = df.iloc[: , 1:] #Delete the first column named "DepMap_ID"

   # threshod_probability = 0.5
    gene_list = list(df.columns.values)
    percentage_list_threshold = []
    for i in gene_list:
        cur_sele = df.loc[df[i] > threshod_probability]
        percentage_list_threshold.append(cur_sele.shape[0]/df.shape[0])

    result = pd.DataFrame({"Gene":gene_list,"probability":percentage_list_threshold})
    
    return(result)

#Step 4: comparison the essentiality between the mutated group and the WT group;
def MDSLP(mut_gene, tumor_type, sample_info, CRISPR_gene_effect, CCLE_mutation, gene_candidates):

    # Example: result = MDSLP("BRCA2", ['pancancer'],  sample_info,  CRISPR_gene_effect, CCLE_mutation, gene_candidates)

    if tumor_type == ['pancancer']:
        pancancer_cls = (sample_info.loc[~sample_info['primary_disease'].isin(['Non-Cancerous','Unknown','Engineered','Immortalized'])])
        pancancer_cls = pancancer_cls.loc[~(pancancer_cls['primary_disease'].isna())]
        
        cl_sele = list(pancancer_cls['DepMap_ID'].values)

    else:
        tumor_selected = tumor_type
        cl_sele = sample_info.loc[sample_info['primary_disease'].isin((tumor_selected))]['DepMap_ID']  
        cl_sele = list(set(list(CRISPR_gene_effect.index.values)).intersection(set(cl_sele)))

    samples_with_mut = set(CCLE_mutation['DepMap_ID'])

    selected_variants = ['Splice_Site',
                         'Frame_Shift_Del',
                         'Frame_Shift_Ins',
                         'Nonstop_Mutation',
                         'In_Frame_Del',
                         'In_Frame_Ins',
                         'Missense_Mutation',
                         'Nonsense_Mutation',
                         'Nonstop_Mutation',
                         'Start_Codon_Del',
                         'Start_Codon_Ins',
                         'Start_Codon_SNP',
                         'Stop_Codon_Del',
                         'Stop_Codon_Del',
                         'Stop_Codon_Ins',
                         'De_novo_Start_OutOfFrame']

    samples_depmap_newname = []
    datatype = "Crispr"
    if datatype == "Crispr":
        #print(CRISPR_gene_effect.head(5))
        samples_depmap = set(CRISPR_gene_effect['DepMap_ID'])
        for sample in samples_depmap:
            samples_depmap_newname.append(sample)

    Samples_with_mut_kd = samples_with_mut.intersection(cl_sele).intersection(samples_depmap_newname)
    CRISPR_gene_effect.index = CRISPR_gene_effect['DepMap_ID']
    CRISPR_gene_effect = CRISPR_gene_effect.iloc[:,1:]
    #print(CRISPR_gene_effect.head(5))
    
    col_names = []
    df = CRISPR_gene_effect
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
    
    
    Depmap_matrix_sele = df.loc[list(Samples_with_mut_kd),gene_candidates].transpose()
    
    
    #geneset = list(gene_candidates)
    #Depmap_matrix_sele = df.loc[:,geneset]
    
    Mut_mat_sele1 = CCLE_mutation.loc[CCLE_mutation['DepMap_ID'].isin(Samples_with_mut_kd)]
    Mut_mat_sele2 = Mut_mat_sele1.loc[Mut_mat_sele1['Variant_Classification'].isin(selected_variants)]
    
    Mut_mat_sele3 = Mut_mat_sele2.loc[Mut_mat_sele2['Hugo_Symbol'] == mut_gene,['Hugo_Symbol','DepMap_ID']]
    
    
    def Cohen_dist(x,y):

        n1 = len(x)
        n2 = len(y)
        s = np.sqrt(((n1 - 1)*(np.std(x))*(np.std(x)) + (n2 - 1) * (np.std(y)) * (np.std(y))) / (n1 + n2 -2))
        d = (np.mean(x) - np.mean(y)) / s
        return(d)
    
    Gene_mut_list = []
    Gene_kd_list = []
    p_list = []
    es_list = []
    size_mut = []
    FDR_List = []
    result = pd.DataFrame()

    p_list_curr = []
    Gene = mut_gene
    Mut_group = list(Mut_mat_sele3.loc[Mut_mat_sele3['Hugo_Symbol'] == Gene]['DepMap_ID'].values)
    #Mut_group = list(Mut_mat_sele3.loc[Mut_mat_sele3['Hugo_Symbol'].isin(mut_gene) ]['DepMap_ID'].values)
    WT_group = list(set(Samples_with_mut_kd) - set(Mut_group))
    print("Number of samples with mutation: " + str(len(Mut_group)))
    print("Number of samples without mutation: " + str(len(WT_group)))
    
    if len(Mut_group) < 5:
        return(pd.DataFrame())
    else:
        for Gene_kd in list(Depmap_matrix_sele.index.values):
            D_mut_new = Depmap_matrix_sele.loc[Gene_kd,Mut_group].values
            D_wt_new = Depmap_matrix_sele.loc[Gene_kd,WT_group].values

            nan_array = np.isnan(D_mut_new)
            not_nan_array = ~ nan_array
            D_mut_new = D_mut_new[not_nan_array]

            nan_array = np.isnan(D_wt_new)
            not_nan_array = ~ nan_array
            D_wt_new = D_wt_new[not_nan_array]

            if len(D_mut_new) > 5:


                Sci_test = stats.ranksums(D_mut_new, D_wt_new, alternative = "less")
                pvalue = Sci_test[1]

                if np.isnan(pvalue) == False:
                    size_mut.append(len(D_mut_new))
                    p_list_curr.append(pvalue)
                    Size_effect =Cohen_dist(D_mut_new, D_wt_new)
                    es_list.append(Size_effect)
                    Gene_mut_list.append(Gene)
                    Gene_kd_list.append(Gene_kd)
        if len(p_list_curr) > 0:

            FDR_List_table = multi.multipletests(p_list_curr, alpha=0.05, method='fdr_bh', is_sorted=False)[1]
            p_list = p_list + p_list_curr
            FDR_List = FDR_List + list(FDR_List_table)

        result = pd.DataFrame({#"Gene_mut": Gene_mut_list, 
                                   #"Gene_mut_symbol": Gene_mut_list,
                                   "Gene_kd": Gene_kd_list, 
                                   #"Gene_kd_symbol":Gene_kd_list,
                                   #"Mutated_samples":size_mut,
                                   "pvalue": p_list_curr, 
                                   "ES":es_list, 
                                   "FDR_by_gene": FDR_List,
                                 #  "FDR_all_exp":FDR_List,
                                 #  "Tumor_type":[','.join(tumor_type)]*len(FDR_List)
                              })
        x = result.loc[result['FDR_by_gene']<0.05]
        x = x.loc[x['ES'] < -0.4]
        x = x.sort_values(by = ['FDR_by_gene'])
        return(x)

#Integrative analysis between mutation dependent essentiality and mutation dependent drug sensitivity
def get_SL_with_multiple_evidence(mut_gene, tumor_type, sample_info, CRISPR_gene_effect, CCLE_mutation, gene_candidates):
    result = MDSLP(mut_gene, tumor_type, sample_info, CRISPR_gene_effect, CCLE_mutation, gene_candidates)
    if result.shape[0] > 1:
        sen_drugs = kg.drug_response_KG_by_gene(mut_gene,"biolink:associated with sensitivity to")
        result_with_both_support = list(set(DT.loc[DT['Drug'].isin(list(sen_drugs['Object_name']))]['Target']).intersection(set(result['Gene_kd'])))
    
        return(result_with_both_support)
    else:
        return([])


