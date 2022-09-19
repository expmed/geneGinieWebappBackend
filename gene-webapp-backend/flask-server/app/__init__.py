from flask import Flask
import pandas as pd
import pdb

app = Flask(__name__)

def create_app():
    """Initialize the flask app instance
    """
    app = Flask(__name__)
    return app

#Read databases and initialize on first request
def read_databases():
        #Read databases
        print("Reading Databases")
        print("Reading gene2pubmed")
        gene2pubmed = pd.read_csv('data/gene2pubmed_primary_genes_with_citations.csv', low_memory=False)
        gene2pubmed.drop(columns=['Unnamed: 0'],inplace=True)
        #Some genes are discontinued and dont have symbols so I remove na
        gene2pubmed = gene2pubmed.dropna()
        #Read gene data for gene cards
        print("Reading gene_data Info")
        # gene_data = pd.read_csv('data/gene_info', delimiter = "\t", low_memory=False)
        #orthologs
        print("Reading Orthologs")
        gene_orthologs = pd.read_csv('data/gene_orthologs', delimiter = "\t")
        #pubmed data paper data (only gene2pubmed ids)
        print("Reading gene2pubmed_papers")
        gene2pubmed_papers = pd.read_csv('data/pubmed_paper_data_gene2pubmed_simplified.csv', low_memory=False)
        gene2pubmed_papers.drop(columns=['Unnamed: 0'],inplace=True)
        print("Finished Reading Databases")
        # pdb.set_trace()
        #Merge databases and processing step
        print("merging databases")
        #Merge both tables I do it on left to keep the rows that I have in gene2pubmed
        gene2pubmed_and_pubmed_paper_data_merged = pd.merge(gene2pubmed,gene2pubmed_papers,left_on='PubMed_ID',right_on='pmid',how='left')
        print(" Finished merging databases")
        #delete pmid column and keep PubMed_ID
        gene2pubmed_and_pubmed_paper_data_merged = gene2pubmed_and_pubmed_paper_data_merged.drop(columns=['pmid'])
        #Remove NA values in publication dates
        gene2pubmed_and_pubmed_paper_data_merged = gene2pubmed_and_pubmed_paper_data_merged.dropna(subset=['pubdate'])
        #setup pubdate as date or int
        gene2pubmed_and_pubmed_paper_data_merged['pubdate'] = gene2pubmed_and_pubmed_paper_data_merged['pubdate'].astype('int')
        print(" Finished processing databases")
        return [gene2pubmed,gene_orthologs,gene2pubmed_papers,gene2pubmed_and_pubmed_paper_data_merged]

#This function is used to add or remove mesh terms for the purposes of the project
#See https://meshb-prev.nlm.nih.gov/treeView go to publication characteristics [V] 
#and then publication formats [v1] for information on mesh terms
# def mesh_terms_mapping:



#Routes Definitions (This route is only used for testing purposes now)
@app.route('/')
def hellow_world():
    return 'flask docker'