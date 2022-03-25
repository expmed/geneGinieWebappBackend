from flask import Flask
import pandas as pd

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
        gene2pubmed = pd.read_csv('data/gene2pubmed_primary_genes.csv', low_memory=False)
        print("Reading Databases")
        gene2pubmed.drop(columns=['Unnamed: 0'],inplace=True)
        print("Reading gene2pubmed")
        #orthologs
        gene_orthologs = pd.read_csv('data/gene_orthologs', delimiter = "\t")
        print("Read orthologs")
        #pubmed data paper data (only gene2pubmed ids)
        gene2pubmed_papers = pd.read_csv('data/pubmed_paper_data_gene2pubmed_simplified.csv', low_memory=False)
        gene2pubmed_papers.drop(columns=['Unnamed: 0'],inplace=True)
        print("Read gene2pubmed_papers")
        print("Finished Reading Databases")
        return [gene2pubmed,gene_orthologs,gene2pubmed_papers]


#Routes Definitions (This route is only used for testing purposes now)
@app.route('/')
def hellow_world():
    return 'flask docker'