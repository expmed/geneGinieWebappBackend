from app import app, read_databases
from flask_sock import Sock
import pdb
import json

import numpy as np
import pandas as pd

from datetime import date
from dateutil.relativedelta import relativedelta

#BioPython
from Bio import Entrez

# app=create_app()

#Routes Definitions (This route is only used for testing purposes now)
@app.route('/')
def test():
    return 'flask docker'

gene2pubmed,gene_orthologs,gene2pubmed_papers, gene2pubmed_and_pubmed_paper_data_merged = read_databases()
#Javascript doesnt like columns named with # so I change the name of the tax_id Column
gene2pubmed_and_pubmed_paper_data_merged.rename(columns={"#tax_id": "tax_id"},inplace=True)

#First we need to filter table by unique 'PubMed_ID' to get a line per publication
unique_gene2pubmed_and_pubmed_paper_data_merged = gene2pubmed_and_pubmed_paper_data_merged.drop_duplicates(subset=['PubMed_ID','GeneID'])
sock = Sock(app)

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
        ### Summary Component data ###
        if data_recieved['type']=='summary_component':
            geneID = data_recieved['msg']
            print('Summary Component Data')
            histogram_all, histogram_all_cum, n_papers, histogram_research, histogram_research_cum, n_research_papers, histogram_reviews, histogram_reviews_cum, n_review_papers = generate_summary_component_data(geneID,unique_gene2pubmed_and_pubmed_paper_data_merged)

            ## Send Summary data ##
            #Send Total Number of papers in all gene2pubmed (number)
            send_data(ws,'summary_component','n_papers_summary_all_gene2pubmed',n_papers)
            #Send Histogram for all gene2pubmed data 
            send_data(ws,'summary_component','histogram_summary_all_gene2pubmed',histogram_all.to_json(orient="records"))
            #Send Cummulative Histogram for only research papers in all gene2pubmed data
            send_data(ws,'summary_component','histogram_cum_summary_all_gene2pubmed',histogram_all_cum.to_json(orient="records"))
            #Send Total Number of papers in all gene2pubmed (number)
            send_data(ws,'summary_component', 'n_papers_summary_research_gene2pubmed',n_research_papers)
            #Send Histogram for only research papers in all gene2pubmed data
            send_data(ws, 'summary_component', 'histogram_summary_research_gene2pubmed',histogram_research.to_json(orient="records"))
            #Send Cummulative Histogram for only research papers in all gene2pubmed data
            send_data(ws, 'summary_component', 'histogram_cum_summary_research_gene2pubmed',histogram_research_cum.to_json(orient="records"))
            send_data(ws,'summary_component', 'n_papers_summary_reviews_gene2pubmed',n_review_papers)
            #Send Histogram for reviews and others papers in all gene2pubmed data
            send_data(ws, 'summary_component', 'histogram_summary_reviews_gene2pubmed',histogram_reviews.to_json(orient="records"))
            #Send Cummulative Histogram for only research papers in all gene2pubmed data
            send_data(ws, 'summary_component', 'histogram_cum_summary_reviews_gene2pubmed',histogram_reviews_cum.to_json(orient="records"))
        ### End Summary Component data ###
        ### Subsection Component data ###
        elif data_recieved['type'] == "subsection_component":
            #Calculate gene counts for all time (Using Primary Genes as Key)
            gene_counts_all_time = bar_chart_counts(gene2pubmed[['primary_gene','Symbol']], 50, 'primary_gene')
            # pdb.set_trace()
            #Send gene counts for all time
            send_data(ws,'subsection_component','gene_counts_all_time',gene_counts_all_time.to_json(orient="records"))
            #Calculate gene counts for last 10 YEARS (Using Primary Genes as Key)
            start_year = time_threshold(10)
            gene_counts_past_10 = unique_gene2pubmed_and_pubmed_paper_data_merged[unique_gene2pubmed_and_pubmed_paper_data_merged['pubdate']>=start_year]
            gene_counts_past_10 = bar_chart_counts(gene_counts_past_10[['primary_gene','Symbol']], 50, 'primary_gene')
            #Send gene counts past 10 years
            send_data(ws,'subsection_component','gene_counts_past_10',gene_counts_past_10.to_json(orient="records"))
            #Calculate gene counts for last 5 YEARS (Using Primary Genes as Key)
            start_year = time_threshold(5)
            gene_counts_past_5 = unique_gene2pubmed_and_pubmed_paper_data_merged[unique_gene2pubmed_and_pubmed_paper_data_merged['pubdate']>=start_year]
            gene_counts_past_5 = bar_chart_counts(gene_counts_past_5[['primary_gene','Symbol']], 50, 'primary_gene')
            #Send gene counts past 5 years
            send_data(ws,'subsection_component','gene_counts_past_5',gene_counts_past_5.to_json(orient="records"))
            #Calculate gene counts for last 1 YEAR (Using Primary Genes as Key)
            start_year = time_threshold(1)
            gene_counts_past_1 = unique_gene2pubmed_and_pubmed_paper_data_merged[unique_gene2pubmed_and_pubmed_paper_data_merged['pubdate']>=start_year]
            gene_counts_past_1 = bar_chart_counts(gene_counts_past_1[['primary_gene','Symbol']], 50, 'primary_gene')
            #Send gene counts past 5 years
            send_data(ws,'subsection_component','gene_counts_past_1',gene_counts_past_1.to_json(orient="records"))
        ### END Subsection Component data ###
        #PAPER DATA SECTION
        elif data_recieved['type'] == "list_component":
            #Store geneID
            geneID = data_recieved['msg']
            geneID = int(geneID)
            #Filter papers with ID
            db = unique_gene2pubmed_and_pubmed_paper_data_merged[unique_gene2pubmed_and_pubmed_paper_data_merged['primary_gene']==geneID]
            #Remove Columns I dont want (MAYBE THIS CHANGES LATER)
            db.pop('tax_id')
            db.pop('authors')
            # db.pop('PubMed_ID')
            db.pop('GeneID')
            db.pop('primary_gene')
            db.pop('Symbol')
            db.pop('mesh_terms')
            db.pop('chemical_list')
            # pdb.set_trace()
            #Order by citations
            db.sort_values(by=['citations'], ascending=False,inplace=True)
            #Past 10
            start_year = time_threshold(10)
            past_10 = db[db['pubdate']>=start_year]
            #Past 5
            start_year = time_threshold(5)
            past_5 = db[db['pubdate']>=start_year]
            #Last Year
            start_year = time_threshold(1)
            past_1 = db[db['pubdate']>=start_year]
            # pdb.set_trace()
            send_data(ws,'subsection_component','all_time',db.to_json(orient="records"))
            send_data(ws,'subsection_component','past_10',past_10.to_json(orient="records"))
            send_data(ws,'subsection_component','past_5',past_5.to_json(orient="records"))
            send_data(ws,'subsection_component','past_1',past_1.to_json(orient="records"))
        #END PAPER DATA SECTION
        
        ### GENE Component data ###
        elif data_recieved['type']=='get_gene_info':
            ### GENE  CARD SECTION ####
            #Store geneID
            geneID = data_recieved['msg']
            #Get gene info data and create gene_object for client (for gene Card)
            print('fetching data')
            gene_info_data = fetch_details([geneID])
            print('creating_object')
            # pdb.set_trace()
            gene_info = create_gene_info_object(gene_info_data[0]) #We should only get one element passing only one ID
            print('gene_info_created')
            #Send gene info data
            send_data(ws,'gene_component', 'gene_info',json.dumps(gene_info))
            #get ranking table and send to client (for gene Card)
            #Calculate gene counts for all time (WE NEED TO POSITION THE TABLE AT THE SPECIFIC GENE)
            # pdb.set_trace()
            gene_counts_all_time = bar_chart_counts(gene2pubmed[['primary_gene','Symbol']], 50, 'primary_gene')
            #Send gene counts for all time
            send_data(ws, 'gene_component', 'gene_counts_all_time',gene_counts_all_time.to_json(orient="records"))
            ### END GENE  CARD SECTION ####


#Function to generate summary_component_data
def generate_summary_component_data(geneID, paper_db):
    """
        db = all database
        research_db = db filtered by research papers
        reviews_db = db filtered by reviews papers
        If we have a geneID we filter db for that gene
        if not we have all genes"""
    if(geneID):
        geneID = int(geneID)
        #Get Paper Data for specific Gene
        db = paper_db[paper_db['primary_gene']==geneID]
    else:
        db = paper_db
    #Filter Research and Review Papers
    research_db = db[db['publication_types'].str.contains('Research')]
    review_db = db[db['publication_types'].str.contains('Review')]
    # pdb.set_trace()
    #Number of papers
    n_papers = len(db)
    #Calculate value counts for histogram all paper types
    histogram_all, histogram_all_cum = calculate_histogram(db['pubdate'])
    #Number of Research papers
    n_research_papers = len(research_db)
    #Calculate value counts for histogram research papers 
    histogram_research, histogram_research_cum = calculate_histogram(research_db['pubdate'])
    #Number of Review papers
    n_review_papers = len(review_db)
    #Calculate  value counts for histogram reviews papers
    histogram_reviews, histogram_reviews_cum = calculate_histogram(review_db['pubdate'])
    # pdb.set_trace()
    return histogram_all, histogram_all_cum, n_papers, histogram_research, histogram_research_cum, n_research_papers, histogram_reviews, histogram_reviews_cum, n_review_papers

#Helper functions for searching GENE INFO in Biopython (Example with gene database)
def search(query,db='gene',sort='relevance',retmax='20'):
    """
        Function to search IDs in NCBI databases using BioPython package
        @requires: from Bio import Entrez
        @params:{
            query: string query to perform search
            db:database we are searching on. (default=gene)
            @returns:Object{
                'Count',
                'RetMax': Max return results
                'RetStart',
                'IdList' : List of Ids in results
                'TranslationSet',
                'TranslationStack',
                'QueryTranslation',
                'WarningList': "Errors in query"
            }
        } 
    """
    Entrez.email = 'agonzamart@gmail.com'
    handle = Entrez.esearch(db=db,
                            sort=sort,
                            retmax=retmax,
                            retmode='xml',
                            term=query)
    results = Entrez.read(handle)
    return results

def fetch_details(id_list,db='gene'):
    """
        AT THE MOMENT I HAVENT FOUND THE NEED TO USE THIS ONE. WE WILL USE IT IN THE SEARCH BAR PROBABLY.
        Function to search gene info in NCBI databases from an array of IDs using BioPython
        @requires: from Bio import Entrez
        @params:{
            id_list: list of Ids to search info on
            db:database we are searching on. (default=gene)
            @returns:Array of Objects [] {
                'Entrezgene_track-info': Basic info. Updated created info.
                'Entrezgene_type': Type
                'Entrezgene_gene': Gene information
                'Entrezgene_prot' : Protein Information. Other names
                'Entrezgene_summary': Description in text
                'Entrezgene_location',
                'Entrezgene_gene-source',
                'Entrezgene_locus',
                'Entrezgene_properties',
                'Entrezgene_comments',
                'Entrezgene_unique-keys',
                'Entrezgene_xtra-index-terms',
                'Entrezgene_xtra-properties'
            }
        } 
    """
    ids = ','.join(id_list)
    Entrez.email = 'agonzamart@gmail.com'
    handle = Entrez.efetch(db=db,
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    results2 = Entrez.parse(handle)
    return results

def create_gene_info_object(gene_data):
    """
        Function to create gene_info_object to send to client
        with the parameters needed for gene_info page
        @params{
            gene_data: Data for a single Gene we want to parse.
    """
    #Create processing for what I want
    gene_info = {
        'Official_Symbol':'',
        'Official_Name':'',
        'Synonyms':'',
        'First_Published':'',
        'Type':'',
        'Organism':'',
        'Lineage':'',
        'Summary':'',
    }
    # pdb.set_trace()
    #Official Symbol
    try:
        gene_info['Official_Symbol'] = gene_data['Entrezgene_gene']['Gene-ref']['Gene-ref_locus']
    except KeyError:
        gene_info['Official_Symbol'] = ''
    #Official Full Name
    try:
        gene_info['Official_Name'] = gene_data['Entrezgene_gene']['Gene-ref']['Gene-ref_desc']
    except KeyError:
        gene_info['Official_Name'] = ''
    #Created On
    try:
        year = gene_data['Entrezgene_track-info']['Gene-track']['Gene-track_create-date']['Date']['Date_std']['Date-std']['Date-std_year']
        month = gene_data['Entrezgene_track-info']['Gene-track']['Gene-track_create-date']['Date']['Date_std']['Date-std']['Date-std_month']
        day = gene_data['Entrezgene_track-info']['Gene-track']['Gene-track_create-date']['Date']['Date_std']['Date-std']['Date-std_day']
        gene_info['First_Published'] = str(month) + '-' + str(day) + '-' + str(year)
    except KeyError:
        gene_info['First_Published'] = ''
    #Also Known as
    try:
        gene_info['Synonyms'] = gene_data['Entrezgene_gene']['Gene-ref']['Gene-ref_syn']
    except KeyError:
        gene_info['Synonyms'] = ''
    #Gene Type
    try:
        gene_info['Type'] = gene_data['Entrezgene_type'].attributes['value']
    except KeyError:
        gene_info['Type'] = ''
    #Organism
    try:
        gene_info['Organism'] = gene_data['Entrezgene_source']['BioSource']['BioSource_org']['Org-ref']['Org-ref_taxname']
    except KeyError:
        gene_info['Organism'] = ''
    #Lineage
    try:
        gene_info['Lineage'] = gene_data['Entrezgene_source']['BioSource']['BioSource_org']['Org-ref']['Org-ref_orgname']['OrgName']['OrgName_lineage']
    except KeyError:
        gene_info['Lineage'] = ''
    #Summary
    try:
        gene_info['Summary'] = gene_data['Entrezgene_summary']
    except KeyError:
        gene_info['Summary'] = ''
    return gene_info

def getCitationCount(publicationIDArray):
    pdb.set_trace()
    publicationIdsStrings = [str(x) for x in publicationIDArray]
    publications_data = fetch_details(publicationIdsStrings)
    citationCounts = [str(publication['PubmedData']['ArticleIdList'][0]) for publication in publications_data['PubmedArticle']]
    return citationCounts

def calculate_histogram(data):
    """
        Function to calculate year histogram data
        Input:data_to_transform to histogram
        Output: Histogram by years, Cumulative Histogram  by years
    """
    MIN_YEAR = 1960
    #calculate counts and order by date
    hist = data.value_counts().sort_index()
    #Start in the MIN DATE
    hist = hist[hist.index>MIN_YEAR]
    # pdb.set_trace()
    hist = hist.to_frame(name="counts")
    #We need to add missing years
    #For that we create an additional DataFrame with all years and merge them
    hist_all_years = pd.DataFrame()
    #range returns sequence -1 so I add a +1 to max
    hist_all_years.index = pd.Series(range(hist.index.min(), hist.index.max()+1))
    hist = hist_all_years.join(hist).fillna(0)
    hist['counts'] = hist['counts'].astype(int)
    #Calculate cumulative histogram
    cum_hist = hist.counts.cumsum().to_frame(name="counts")
    #Add a column for the years
    hist['ID'] = hist.index
    cum_hist['ID'] = cum_hist.index
    return hist,cum_hist

def time_threshold(years):
    """
        Function to calculate date start thresholds for data year filters
        Input: num_years (int)
        Return: year (int)
    """
    relative_date = date.today() - relativedelta(years=years)
    threshold_year_start = int(relative_date.strftime('%Y'))
    return threshold_year_start


def bar_chart_counts(data, num_results,index_col):
    """
        Function that calculates counts for a column
        Inputs:
            - data: Dataframe and column to calculates counts from
            - num_results to return from function. Returns top num_results
        Output: data_frame with column counts
    """
    column_order = ['Ranking',index_col,'Symbol','counts']
    counts_dataFrame = data[index_col].value_counts().to_frame(name='counts')
    #Reset index
    counts_dataFrame = counts_dataFrame.reset_index()
    #Rename column for the IDs
    counts_dataFrame = counts_dataFrame.rename(columns={"index": index_col})
    #Merge with Symbols
    table_with_symbols_with_duplicates = pd.merge(counts_dataFrame,data,on=index_col,how='left')
    #Clean merged table of duplicates (select last to get capitalized Version of Symbol) and reset index
    counts_dataFrame = table_with_symbols_with_duplicates.drop_duplicates(subset=[index_col], keep='last')
    #Order by counts and store the top 50 to send
    counts_dataFrame = counts_dataFrame.sort_values(by=['counts'], ascending=False)[:num_results]
    #Remove current index to create ranking column new index is the ranking
    #First reset the current index
    counts_dataFrame.reset_index(drop=True,inplace=True)
    #then reset again to store the ranking
    counts_dataFrame.reset_index(inplace=True)
    counts_dataFrame.rename(columns={"index": "Ranking"},inplace=True)
    #We add 1 to start at 1 instead of at 0.
    counts_dataFrame['Ranking'] = counts_dataFrame['Ranking'] + 1
    #Reorder columns
    counts_dataFrame = counts_dataFrame[column_order]
    # pdb.set_trace()
    return counts_dataFrame

def send_data(socket,component,type,data):
    """ Function to send data back to client """
    print('===============')
    print(component)
    print(type)
    print('sending data')
    print('===============')
    msg = {
        'component':component,
        'type':type,
        'data':data
    }
    socket.send(json.dumps(msg))
    