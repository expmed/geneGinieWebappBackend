from app import create_app
from flask_sock import Sock
import pdb
import json

import numpy as np
import pandas as pd

app=create_app()

sock = Sock(app)

#Read databases
print("Reading Databases")
gene2pubmed = pd.read_csv('data/gene2pubmed_primary_genes.csv', low_memory=False)
gene2pubmed.drop(columns=['Unnamed: 0'],inplace=True)
#orthologs
gene_orthologs = pd.read_csv('data/gene_orthologs', delimiter = "\t")

#pubmed data paper data (only gene2pubmed ids)
gene2pubmed_papers = pd.read_csv('data/pubmed_paper_data_gene2pubmed_simplified.csv', low_memory=False)
gene2pubmed_papers.drop(columns=['Unnamed: 0'],inplace=True)
print("Finished Reading Databases")

# @app.route("/", methods=['get'])
# def home():
#     # print (gene2pubmed)
#     return "Hello world"

#Socket definition
@sock.route('/echo')
def echo(ws):
    while True:
        data_recieved = 0
        data_recieved = ws.receive()
        data_recieved = json.loads(data_recieved)
        # pdb.set_trace()
        if data_recieved['type']=='overview':
            print(data_recieved)
            send_data(ws,'overview',gene2pubmed_papers[:10].to_json(orient="records"))
        elif data_recieved['type']=='search_gene':
            print('search_gene')
            pubmed_gene_publications = search_gene(int(data_recieved['msg']))
            #send table data
            send_data(ws,'table_data',pubmed_gene_publications.to_json(orient='records'))
            #Calculate value counts for histogram
            # pdb.set_trace()
            histogram_data = calculate_date_histogram(pubmed_gene_publications['pubdate'])
            # pdb.set_trace()
            #send hist data
            send_data(ws,'hist_data', histogram_data.to_json(orient='records'))

"""Function that finds papers in PubMed for a specific gene passing a geneID"""
def search_gene(geneID):
    #get primary gene
    primary_gene = gene2pubmed[gene2pubmed['GeneID']==geneID]['primary_gene'].values[0]
    #Get orthologs list
    list_orthologs = gene_orthologs[gene_orthologs['GeneID']==primary_gene]['Other_tax_id'].values
    #append gene user searched for
    list_orthologs = np.append(list_orthologs,primary_gene)
    #Get all PubMed_ID for these genes
    pubmed_id_list = gene2pubmed[gene2pubmed['GeneID'].isin(list_orthologs)]['PubMed_ID'].unique()
    #Search publications in gene2pubmed_papers
    papers_for_gene = gene2pubmed_papers[gene2pubmed_papers['pmid'].isin(pubmed_id_list)]
    #Drop columns we dont want
    papers_for_gene.pop('mesh_terms')
    papers_for_gene.pop('chemical_list')
    return papers_for_gene

def calculate_date_histogram(data):
    #calculate counts and order by date
    table = data.value_counts().sort_index()
    table = table.to_frame()
    #Add a column for the years
    table['year'] = table.index
    return table


"""Function to send data back to client"""
def send_data(socket,type,data):
    print('sending data')
    msg = {
        'type':type,
        'data':data
    }
    socket.send(json.dumps(msg))
    