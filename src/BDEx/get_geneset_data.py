import pandas as pd
import os

def get_geneset(source, pathway_name):
    genesets = []
    if source == "msigdb":
        url = "https://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=" + pathway_name + "&fileType=txt"
        genesets = pd.read_csv(url)
        genesets = list(genesets[pathway_name].values)[1:len(genesets)]
    else:
        print("Usage: get_geneset('msigdb', pathway_name)")

    return(genesets)

def get_data(Dataset, tableName):
    temp_dir = "./Temp/"
    path_to_file = temp_dir + tableName + '.csv'
    if os.path.exists(temp_dir) == False:
        try:
            os.makedirs(temp_dir)
        except OSError:
            print ("Creation of the directory %s failed" % temp_dir)
        else:
            print ("Successfully created the directory %s " % temp_dir)
    else:
        print ("INfO:  %s already exists!" % temp_dir)
    
    from os.path import exists
    if os.path.exists(path_to_file):
        result = pd.read_csv(path_to_file)
        print("Read file from " + path_to_file)
        return(result)
    else:
        resource_table = pd.read_csv("../../data/resource_table.csv")
        #print(resource_table['Dataset'])
        url = ''
        if Dataset in set(resource_table['Dataset'].values):
            sele = resource_table.loc[resource_table['Dataset'] == Dataset]
            if tableName in set(sele['Table'].values):
                sele1= sele.loc[sele['Table'] == tableName]
                url = sele1['URL'].values[0]
                print(url)
            else:
                print("The following tables are provided:")
                for table in sele['Table']:
                    print(table)
        else:
            print("The following Datasets are supported:")
            for i in set(resource_table['Dataset']):
                print(i)

        if len(url) > 0:
            result = pd.read_csv(url)
            print("Read file from remote resource.")
            result.to_csv(path_to_file, index=False)
            return(result)
        else:
            return()
