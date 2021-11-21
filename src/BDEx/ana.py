from scipy import stats
import pandas as pd
import numpy as np

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