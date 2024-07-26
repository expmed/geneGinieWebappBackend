from flask import Flask
from flask_cors import CORS
from sqlalchemy.sql import text
import sqlalchemy as db
import pandas as pd
import pdb

app = Flask(__name__)
#Allow CORS in the APP (I think this was needed for the database to work)
CORS(app, origins=['http://localhost:4200'])

# Define the path to the SQLite database file
db_path = 'data/ginie.db' 
engine = db.create_engine(f'sqlite:///{db_path}')

def create_app():
    """Initialize the flask app instance
    """
    app = Flask(__name__)
    return app


"""Helper function to get values from the database SQLITE3 NOW"""
def get_data(db_name,query=''):
    #No query: Get all database
    if query =='':
        connection = engine.connect()
        result = pd.read_sql_table(db_name,con=connection)
        connection.close()
    else:
        sql = f'{query}'
        with engine.connect().execution_options(autocommit=True) as conn:
            query = conn.execute(text(sql))         
        result = pd.DataFrame(query.fetchall())
    return result
    

#Read databases and initialize on first request
def read_databases():
        print('Reading Databases')
        # gene2pubmed_db = pd.read_feather('data/gene2pubmed_p_genes_db.feather')
        gene2pubmed_db = get_data('gene2pubmed')
        gene2pubmed_papers_db =  get_data('gene2pubmed_papers')
        # gene2pubmed_papers_db =  pd.read_csv('data/gene2pubmed_pub_db.csv')
        authors_db =  get_data('authors')
        # gene_info_db = pd.read_feather('data/gene_info.feather')
        gene_info_db = get_data('gene_info')
        symbol_group = get_data('symbol_group')
        #transform to proper lists
        symbol_group['idList']= symbol_group['idList'].apply(lambda x: eval(x))
        #transform to list of ints
        symbol_group['idList'] = symbol_group['idList'].apply(lambda x: [int(i) for i in x])
        citations_db = pd.read_feather('data/citations_db.feather')
        print('Finished Reading Databases')
        print('start Processing for databases')
        
        #Format authors
        # gene2pubmed_papers_db['authors'].fillna("",inplace=True)
        # gene2pubmed_papers_db['authors'] = gene2pubmed_papers_db['authors'].apply(lambda x: x.replace('|','_')).apply(lambda x: x.split(';'))
        # #Fill synonyms empty values and rename #tax_id column to tax_id
        # gene_info_db['Synonyms'] = gene_info_db['Synonyms'].fillna('')
        # gene_info_db.rename(columns={"#tax_id": 'tax_id'},inplace=True)

        #Filter gene Info database for only primary genes (for now at least I dont see reason to have everything)
        # primary_genes_IDs = gene2pubmed_db['primary_gene'].unique()
        # gene_info_db = gene_info_db[gene_info_db['GeneID'].isin(primary_genes_IDs)]
        # primary_genes_symbols = gene_info_db['Symbol'].unique()

        # print('End Processing for databases')

        # print('genes home subsection')

        # gene_home = pd.merge(gene2pubmed_db,gene_info_db[['GeneID','Symbol']],left_on='primary_gene',right_on='GeneID',how='left').drop_duplicates()
        # #Attach date
        # gene_home = pd.merge(gene_home,gene2pubmed_papers_db[['PubMed_ID','pubdate']],left_on='PubMed_ID',right_on='PubMed_ID',how='left').drop_duplicates()
        # gene_home.drop(columns=['tax_id','GeneID_x','GeneID_y'],inplace=True)
        # gene_home = gene_home.drop_duplicates()
        # #Drop na values (MIGHT NEED TO REVIEW LATER - Some symbols and dates are nan)
        # gene_home = gene_home.dropna()
        # #set date column as int
        # gene_home['pubdate'] = gene_home['pubdate'].astype('int')

        # print('end genes home subection')

        return [gene2pubmed_db,gene2pubmed_papers_db,gene_info_db,authors_db,citations_db,symbol_group]



#Routes Definitions (This route is only used for testing purposes now)
@app.route('/')
def hellow_world():
    return 'flask docker'