import pandas as pd
import requests
import json



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
