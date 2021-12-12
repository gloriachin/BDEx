from scipy import stats
import pandas as pd
import numpy as np
from scipy.stats import distributions
import statsmodels.stats
from statsmodels.stats import multitest

def rank_simple(vector):
    return sorted(range(len(vector)), key=vector.__getitem__)


def rankdata(a):
    n = len(a)
    ivec=rank_simple(a)
    svec=[a[rank] for rank in ivec]
    sumranks = 0
    dupcount = 0
    newarray = [0]*n
    for i in range(0,n):
        sumranks += i
        dupcount += 1
        if i==n-1 or svec[i] != svec[i+1]:
            averank = sumranks / float(dupcount) + 1
            for j in range(i-dupcount+1,i+1):
                newarray[ivec[j]] = averank
            sumranks = 0
            dupcount = 0
    return newarray

def rp_spearman_batch(x_ranks, y_ranks):
    dist_array = x_ranks - y_ranks
    d2 = np.einsum("ij,ij->i", dist_array, dist_array, optimize=True)
    n = x_ranks.shape[1]
    r = 1 - 6 * d2 / (n*(n*n-1))
    fac = (r + 1) * (1 - r)
    t = r * np.sqrt((n-2)/fac)
    t = np.where(fac <= 0.0, 0.0, t)
    return (r, n, t)


def measure_cor_spearman_batch(dic_expr_rank, genes_for_query1, genes_for_query2):
    gene1_list = []
    gene2_list = []
    r_list = []
    p_list = []

    g1_list = []
    g2_list = []

    def batch_calculate():
        batch_size = len(g1_list)
        if batch_size == 0:
            return
        rank1 = [dic_expr_rank[g] for g in g1_list]
        rank2 = [dic_expr_rank[g] for g in g2_list]
        r_arr, n, t_arr = rp_spearman_batch(np.array(rank1), np.array(rank2))
    
        p_arr = stats.t.sf(np.abs(t_arr), df=n-2)*2
        r_idx = np.where(np.abs(r_arr) > 0, True, False)
        p_idx = np.where(p_arr < 0.05, True, False)
        idx = np.logical_and(r_idx, p_idx)
        idx = np.where(idx)
        r_list.extend(np.round(r_arr[idx], 3).tolist())
        p_list.extend(p_arr[idx].tolist())
        gene1_list.extend(np.array(g1_list)[idx].tolist())
        gene2_list.extend(np.array(g2_list)[idx].tolist())
        g1_list.clear()
        g2_list.clear()

    for gene1 in genes_for_query1:
        for gene2 in genes_for_query2:
            if gene1 < gene2:
                g1_list.append(gene1)
                g2_list.append(gene2)
                if len(g1_list) >= 50000:
                    batch_calculate()

    batch_calculate()

    result = pd.DataFrame({"Gene1": gene1_list,
                            "Gene2": gene2_list,
                            "rho":r_list,
                            "pvalue":p_list,
                            })

    return(result)


def pd2dic(df):
    col_names = []
    if "Unnamed: 0" in df.columns:
        for col in df.columns:
            if col != "Unnamed: 0":
                col_names.append(col)

        df.drop(columns=df.columns[0], 
                axis=1, 
                inplace=True)
    else:
        for col in df.columns:
            if col != "Unnamed: 0":
                col_names.append(col)

    df.columns = col_names

    result_dic = {}
    for col in list(df.columns):
        a = df[col].values
        rank_a = rankdata(a)
        result_dic[col] = rank_a


    return(result_dic)


def check_gene_groups_from_correlation(correlation_df,threshold_cor):
    group1 = set() #postively correlated group
    group2 = set() #negatively correlated group
    group3 = set() #ambiguious correlated group

        #step1:check positive correlation:
    for i in range(0, correlation_df.shape[0]):
        if correlation_df.iloc[i,2] >threshold_cor and correlation_df.iloc[i,3] < 0.05:
            if len(group1) == 0:
                group1.add(correlation_df.iloc[i,0])
                group1.add(correlation_df.iloc[i,1])
            else:
                if correlation_df.iloc[i,0] in group1:
                    group1.add(correlation_df.iloc[i,1])
                if correlation_df.iloc[i,1] in group1:
                    group1.add(correlation_df.iloc[i,0])

        ##step2:check negative correlation:
    for i in range(0, correlation_df.shape[0]):
        if correlation_df.iloc[i,2] < -1 * threshold_cor and correlation_df.iloc[i,3] < 0.05:
            if len(group1) == 0:
                group1.add(correlation_df.iloc[i,0])
                group2.add(correlation_df.iloc[i,1])
            else:
                if correlation_df.iloc[i,0] in group1:
                    group2.add(correlation_df.iloc[i,1])
                if correlation_df.iloc[i,1] in group1:
                    group2.add(correlation_df.iloc[i,0])


        ##step3:check not significant associated correlation:
    for i in range(0, correlation_df.shape[0]):
        if np.abs(correlation_df.iloc[i,2]) < threshold_cor:
            if len(group1) == 0 and len(group2) == 0:
                print("No significantly correlated genes with thereshold")

            else:
                if correlation_df.iloc[i,0] in group1:
                    group3.add(correlation_df.iloc[i,1])
                elif correlation_df.iloc[i,1] in group1:
                    group3.add(correlation_df.iloc[i,0])
                elif correlation_df.iloc[i,1] in group2:
                    group3.add(correlation_df.iloc[i,0])
                elif correlation_df.iloc[i,0] in group2:
                    group3.add(correlation_df.iloc[i,1])
                else:
                    group3.add(correlation_df.iloc[i,0])
                    group3.add(correlation_df.iloc[i,1])
    group3 = group3 - group1 - group2
    return([list(group1), list(group2), list(group3)])


def _normtest_finish(z, alternative):
    """Common code between all the normality-test functions."""
    if alternative == 'less':
        prob = distributions.norm.cdf(z)
    elif alternative == 'greater':
        prob = distributions.norm.sf(z)
    elif alternative == 'two-sided':
        prob = 2 * distributions.norm.sf(np.abs(z))
    else:
        raise ValueError("alternative must be "
                         "'less', 'greater' or 'two-sided'")
    if z.ndim == 0:
        z = z[()]
    return z, prob

def rank_sum_test_batch_CCLE(mut_samples, sele_expr_ranked, alternative):
    '''
    mut_samples: a selection of samples from the list
    sele_expr_ranked: the expression of genes, two columns are included: 
    1) sample ids; 2)expression values; index are ranked from 1 to N, where 1 means lower expression.
    alternative: 'less' or 'greater' or 'two-sided'
    '''
    x_rank = list(sele_expr_ranked.loc[sele_expr_ranked['Unnamed: 0'].isin(mut_samples)].index)
    s = np.sum(x_rank, axis=0)
    n1 = len(mut_samples)
    n2 = sele_expr_ranked.shape[0] - n1
    expected = n1 * (n1 + n2 +1) / 2.0
    z = (s - expected) / np.sqrt(n1*n2*(n1+n2+1)/12.0)
    z, prob = _normtest_finish(z, alternative)
    return([n1,n2,z,prob])

def mut_dep_expr_CCLE(Sample_info, CCLE_mut, CCLE_expr, tumortype, alternative, mut_list, expr_list):
    selected_disease = [tumortype] #Selection of tumor types to be analized.
    selected_samples = list(set(Sample_info.loc[Sample_info['primary_disease'].isin(selected_disease) ]['DepMap_ID']))
    selected_mut = CCLE_mut.loc[CCLE_mut['DepMap_ID'].isin(selected_samples)]

    CCLE_expr.index = CCLE_expr['Unnamed: 0']
    selected_expr = CCLE_expr.loc[CCLE_expr['Unnamed: 0'].isin(selected_samples)]

    All_samples_with_m_e = list(set(selected_mut['DepMap_ID']).intersection(set(selected_expr['Unnamed: 0']))) #Selection of samples with both mutation data and expression data.

    CCLE_functional = selected_mut.loc[selected_mut['Variant_Classification'].isin(['De_novo_Start_OutOfFrame',
                                                                            'Frame_Shift_Del',
                                                                            'Frame_Shift_Ins',
                                                                            'Missense_Mutation',
                                                                            'Nonsense_Mutation',
                                                                            'Nonstop_Mutation',
                                                                            'Splice_Site',
                                                                            'Start_Codon_Del',
                                                                            'Start_Codon_Ins',
                                                                            'Start_Codon_SNP',
                                                                            'Stop_Codon_Del',
                                                                            'Stop_Codon_Ins'
                                                                         ])] #Selection of potential functional genetical alterations

    CCLE_functional = CCLE_functional.loc[CCLE_functional['DepMap_ID'].isin(All_samples_with_m_e)]
    
    # Start the analysis
    result_expr = []
    result_mut = []
    result_mut_size = []
    result_wt_size = []
    result_z = []
    result_p = []

    dic_mut_samples = {}
    for i in range(0,CCLE_functional.shape[0]):
        symbol = CCLE_functional.iloc[i,0]
        sample = CCLE_functional.iloc[i,15]
        if symbol in dic_mut_samples:
            dic_mut_samples[symbol].add(sample)
        else:
            dic_mut_samples[symbol] = set([sample])
                
    for gene_expr in set(list(set(selected_expr.columns) - set(["Unnamed: 0"]))).intersection(expr_list):
        sele_expr = selected_expr.loc[All_samples_with_m_e,["Unnamed: 0",gene_expr]]
        sele_expr_ranked = sele_expr.sort_values(by = [gene_expr])
        sele_expr_ranked.index = range(1,sele_expr_ranked.shape[0]+1) #we will rank the data dataframe by gene expression first, so we don't need to rank it for multiple times

        for mut in set(dic_mut_samples.keys()).intersection(mut_list):
            if len(dic_mut_samples[mut]) > 3:
                mut_samples = dic_mut_samples[mut]
                result = rank_sum_test_batch_CCLE(mut_samples, sele_expr_ranked, alternative)
                #result[3] < 0.05:
                result_expr.append(gene_expr.split(' ')[0])
                result_mut.append(mut)
                result_mut_size.append(result[0])
                result_wt_size.append(result[1])
                result_z.append(result[2])
                result_p.append(result[3])

    Final = pd.DataFrame({"Gene_expr":result_expr,
                          "Gene_mut":result_mut,
                          "mut_size":result_mut_size,
                          "wt_size":result_wt_size,
                          "z":result_z,
                          "p":result_p
                         })
    if len(list(Final['p'])) > 0:
        fdr_list = multitest.multipletests(list(Final['p']),  alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)[1]
        Final['fdr'] = fdr_list
    Final['Disease'] = [tumortype] * Final.shape[0]
    return(Final)