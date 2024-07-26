from app import app, read_databases
from flask import request, jsonify
from flask_sock import Sock

import pdb
import json

import numpy as np
import pandas as pd

import re

from datetime import date
from dateutil.relativedelta import relativedelta

from urllib.error import HTTPError

#BioPython
from Bio import Entrez

gene2pubmed_db,gene2pubmed_papers_db,gene_info_db,authors_db,citations_db,symbol_group = read_databases()

sock = Sock(app)

#Socket definition
@sock.route('/echo')
def echo(ws):
    while True:
        data_recieved = 0
        data_recieved = ws.receive()
        data_recieved = json.loads(data_recieved)
        #store get parametersinder
        component_type = data_recieved['component_type']
        data_filter = data_recieved['value']
        data_options = data_recieved['options']

        ### Seach Results Component data ###
        if component_type =='search-results-component':
            search_results = []
            print('search results component data')
            ##Filter gene db
            #create search term with up and lower case
            # searchTerm = [data_filter.lower(),data_filter.upper()]
            #Primarty genes
            # symbol_results = gene_info_db[gene_info_db['Symbol'].str.contains('|'.join(searchTerm))]
            #Add synonyms
            # synosyms_results = gene_info_db[gene_info_db['Synonyms'].str.contains('|'.join(searchTerm))]
            #concat results
            # search_results = pd.concat([symbol_results,synosyms_results]).drop_duplicates()
            #Search for term in API
            search_results_from_API = search(data_filter)
            #Get list of IDs
            ID_results = search_results_from_API['IdList']
            #Transform to ints
            ID_results = [int(i) for i in ID_results]
            #Filter db with those IDs
            search_results = gene_info_db[gene_info_db['GeneID'].isin(ID_results)]
            #Order search results based on API Order: 
            #1. GeneID to categories with order specified by ID_results
            search_results['GeneID'] = pd.Categorical(search_results['GeneID'],ID_results)
            #2. Tell to sort by GENEID it applies the order specified in ID_results
            search_results = search_results.sort_values('GeneID')
            #Clean db columns
            search_results.drop(columns=[
                'Nomenclature_status',
                'Modification_date',
                'Feature_type',
                'dbXrefs',
                'map_location',
                'Symbol_from_nomenclature_authority',
                'LocusTag',
                'Full_name_from_nomenclature_authority'
                ],inplace=True)
            #calculate symbol size to order by shortest text
            search_results['Symbol_length'] = search_results['Symbol'].str.len()
            search_results = search_results.sort_values(by=['Symbol_length'])
            search_results.pop('Symbol_length')
            send_data(ws,'search-results-component','genes',search_results.to_json(orient="records"))
        ### Summary Component data ###
        elif component_type == 'summary_component':
            page = data_options[0]
            #Filter db based on type to show only this author for author page
            if page == 'author':
                authorName = data_filter
                print(authorName) 
                author_filtered_db = authors_db[authors_db['author'].str.contains(authorName)]
                pubmedIDs = author_filtered_db['PubMed_ID'].values.tolist()
                if(len(pubmedIDs)==0):
                    key = 'GeneID'
                    pubmedIDs = gene2pubmed_db[gene2pubmed_db[key]==int(data_filter)]['PubMed_ID'].values.tolist()
                db = gene2pubmed_db[gene2pubmed_db['PubMed_ID'].isin(pubmedIDs)]
            else:
                db = gene2pubmed_db
                #If we have a geneID we filter db by gene.
                if(data_filter):
                    geneID = int(data_filter)
                    #Get Paper Data for specific Gene
                    db = db[db['Primary_gene']==geneID]
                    #If we dont find data. Non-primary gene
                    if(len(db)==0):
                            key = 'GeneID'
                            db = gene2pubmed_db[gene2pubmed_db[key]==geneID]
            print('Summary Component Data')
            histogram_all, histogram_all_cum, n_papers, histogram_research, histogram_research_cum, n_research_papers, histogram_reviews, histogram_reviews_cum, n_review_papers = generate_summary_component_data(db)


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
        elif component_type == "subsection_component":
            #We play with the value of key to return primary genes or geneIds
            key = ''
            #Store page and data type for subsection
            data_type = data_options[0]
            page = data_options[1]
            #Select db based on page data type and filters
            db = 0
            paper_db = gene2pubmed_papers_db
            if page == 'gene' or page =='home':
                if data_type == 'gene':
                    key = 'Symbol_group'
                    #We need to remove the duplicated pubids for the counts
                    #because multiple smaller genes can be attached to same pubid all assigned to same symbol_group
                    db = gene2pubmed_db.drop_duplicates(subset=['PubMed_ID','Symbol_group'])
                    if data_filter:
                        data_filter = int(data_filter)
                        #We filter 
                        db = gene2pubmed_db[gene2pubmed_db[key]==data_filter].drop_duplicates(subset=['PubMed_ID','Symbol_group'])
                         #If it is empty it is not a symbol_group so we search on geneID (NOT tested yet since new change) ??
                        if(len(db)==0):
                            key = 'GeneID'
                            db = gene2pubmed_db[gene2pubmed_db[key]==data_filter]
                elif data_type == 'author':
                    key = 'Primary_gene'
                    db = authors_db
                    if data_filter:
                        data_filter = int(data_filter)
                        pubmed_ids = gene2pubmed_db[gene2pubmed_db[key]==data_filter]['PubMed_ID'].values.tolist()
                        if(len(pubmed_ids)==0):
                            key = 'GeneID'
                            pubmed_ids = gene2pubmed_db[gene2pubmed_db[key]==data_filter]['PubMed_ID'].values.tolist()
                        #Get the authors
                        db = authors_db[authors_db['PubMed_ID'].isin(pubmed_ids)]
            elif page == 'author':
                author_filtered_db = authors_db[authors_db['author'].str.contains(data_filter)]
                pubmedIDs = author_filtered_db['PubMed_ID'].values.tolist()
                paper_db = paper_db[paper_db['PubMed_ID'].isin(pubmedIDs)]
                if data_type == 'gene':
                    # db = gene2pubmed_db[gene2pubmed_db['PubMed_ID'].isin(pubmedIDs)]
                    db = gene2pubmed_db[gene2pubmed_db['PubMed_ID'].isin(pubmedIDs)]
                elif data_type == 'author':
                    db = authors_db
                    #We need the filter author table by publications in paper_db published by the author
                    pubmedIDs = paper_db['PubMed_ID'].values.tolist()
                    #search for the authors based on the publication IDs
                    author_filtered_db = db[db['PubMed_ID'].isin(pubmedIDs)]
                    db = author_filtered_db
                    #Remove self author from list           
            #Process and send based on page
            if data_type == 'gene':
                ## GENE PUBLICATION COUNTS ##
                #Calculate gene counts for all time
                gene_counts_all_time = db[[key]]
                gene_counts_all_time = bar_chart_counts(gene_counts_all_time, 50, key)
                #Send gene counts for all time
                send_data(ws,'subsection_component','gene_counts_all_time',gene_counts_all_time.to_json(orient="records"))
                #Calculate gene counts for last 10 YEARS (Using Primary Genes as Key)
                start_year = time_threshold(10)
                gene_counts_past_10 = db[db['Pubdate']>=start_year][[key]]
                gene_counts_past_10 = bar_chart_counts(gene_counts_past_10, 50, key)
                #Send gene counts past 10 years
                send_data(ws,'subsection_component','gene_counts_past_10',gene_counts_past_10.to_json(orient="records"))
                #Calculate gene counts for last 5 YEARS (Using Primary Genes as Key)
                start_year = time_threshold(5)
                gene_counts_past_5 = db[db['Pubdate']>=start_year][[key]]
                gene_counts_past_5 = bar_chart_counts(gene_counts_past_5, 50, key)
                #Send gene counts past 5 years
                send_data(ws,'subsection_component','gene_counts_past_5',gene_counts_past_5.to_json(orient="records"))
                #Calculate gene counts for last 1 YEAR (Using Primary Genes as Key)
                start_year = time_threshold(1)
                gene_counts_past_1 = db[db['Pubdate']>=start_year][[key]]
                gene_counts_past_1 = bar_chart_counts(gene_counts_past_1, 50, key)
                #Send gene counts past 1 years
                send_data(ws,'subsection_component','gene_counts_past_1',gene_counts_past_1.to_json(orient="records"))
                ## END GENE PUBLICATION COUNTS ##
            elif data_type == 'author':
                authorkey = ['LastName','ForeName','Initials'] 
                max_n_elements = 200
                #Calculate author counts for all time
                author_counts_all_time = bar_chart_counts(db,max_n_elements,authorkey)
                #Send author counts for all time
                send_data(ws,'subsection_component','author_counts_all_time',author_counts_all_time.to_json(orient="records"))
                #Calculate author counts for last 10 YEARS
                start_year = time_threshold(10)
                # 1. Get pubmedIds with correct threshold
                pubmed_ids = paper_db[paper_db['Pubdate']>=start_year]['PubMed_ID'].values.tolist()
                # 2. Filter author db with those Ids.
                author_counts_past_10 = db[db['PubMed_ID'].isin(pubmed_ids)]
                # 3. Get the counts
                author_counts_past_10 = bar_chart_counts(author_counts_past_10,max_n_elements,authorkey)
                #Send Author counts past 10 Years
                send_data(ws,'subsection_component','author_counts_past_10',author_counts_past_10.to_json(orient="records"))
                #Calculate author counts for last 5 YEARS
                start_year = time_threshold(5)
                # 1. Get pubmedIds with correct threshold
                pubmed_ids = paper_db[paper_db['Pubdate']>=start_year]['PubMed_ID'].values.tolist()
                # 2. Filter author db with those Ids.
                author_counts_past_5 = db[db['PubMed_ID'].isin(pubmed_ids)]
                author_counts_past_5 = bar_chart_counts(author_counts_past_5,max_n_elements,authorkey)
                #Send Author counts past 5 Years
                send_data(ws,'subsection_component','author_counts_past_5',author_counts_past_5.to_json(orient="records"))
                #Calculate author counts for last 1 YEARS
                start_year = time_threshold(1)
                # 1. Get pubmedIds with correct threshold
                pubmed_ids = paper_db[paper_db['Pubdate']>=start_year]['PubMed_ID'].values.tolist()
                # 2. Filter author db with those Ids.
                author_counts_past_1 = db[db['PubMed_ID'].isin(pubmed_ids)]
                author_counts_past_1 = bar_chart_counts(author_counts_past_1,max_n_elements,authorkey)
                #Send Author counts past 1 Years
                send_data(ws,'subsection_component','author_counts_past_1',author_counts_past_1.to_json(orient="records"))
                ## END AUTHOR PUBLICATION COUNTS ##
        ### END Subsection Component data ###
        #PAPER DATA SECTION
        elif component_type == "list_component":
            page = data_options[0]
            if page == 'gene':
                #Store geneID
                geneID = int(data_recieved['value'])
                #Filter publication Ids with Gene
                pubmed_ids = gene2pubmed_db[gene2pubmed_db['Primary_gene']==geneID]['PubMed_ID'].values.tolist()
                #IF empty we are searching on a non-primary gene
                if(len(pubmed_ids)==0):
                    key = 'GeneID'
                    pubmed_ids = gene2pubmed_db[gene2pubmed_db[key]==int(data_filter)]['PubMed_ID'].values.tolist()
            elif page== 'author':
                #store authorId
                auhorID = data_recieved['value']
                #Filter publication Ids with Author
                pubmed_ids = authors_db[authors_db['author'].str.contains(auhorID)]['PubMed_ID'].values.tolist()
            #Get the papers
            paper_db = gene2pubmed_papers_db[gene2pubmed_papers_db['PubMed_ID'].isin(pubmed_ids)]
            #Remove Columns I dont want (MAYBE THIS CHANGES LATER)
            paper_db.drop(columns=[
                'Publication_types',
                'Mesh_terms',
                'Chemical_list'
                ],inplace=True)
            #Add Citations. NaN citations are assumed to be 0 
            paper_db = pd.merge(paper_db,citations_db,left_on='PubMed_ID',right_on='PubMed_ID',how='left').fillna(0)
            paper_db = paper_db.astype({'citations':'int'})
            #Order by citations
            paper_db.sort_values(by=['citations'], ascending=False,inplace=True)
            #Past 10
            start_year = time_threshold(10)
            past_10 = paper_db[paper_db['Pubdate']>=start_year]
            #Past 5
            start_year = time_threshold(5)
            past_5 = paper_db[paper_db['Pubdate']>=start_year]
            #Last Year
            start_year = time_threshold(1)
            past_1 = paper_db[paper_db['Pubdate']>=start_year]
            send_data(ws,'subsection_component','all_time',paper_db.to_json(orient="records"))
            send_data(ws,'subsection_component','past_10',past_10.to_json(orient="records"))
            send_data(ws,'subsection_component','past_5',past_5.to_json(orient="records"))
            send_data(ws,'subsection_component','past_1',past_1.to_json(orient="records"))
        #END PAPER DATA SECTION
        
        ### GENE Component data ###
        elif component_type =='get_gene_info':
            ### GENE  CARD SECTION ####
            #Store geneID
            geneID = data_recieved['value']
            #Get gene info data and create gene_object for client (for gene Card)
            print('fetching data')
            gene_info_data = fetch_details([geneID])
            print('creating_object')
            gene_info = create_gene_info_object(gene_info_data[0],geneID) #We should only get one element passing only one ID
            print('gene_info_created')
            #Send gene info data
            send_data(ws,'gene_component', 'gene_info',json.dumps(gene_info))
            #get ranking table and send to client (for gene Card)
            #Calculate gene counts for all time (WE NEED TO POSITION THE TABLE AT THE SPECIFIC GENE)
        ### END GENE  CARD SECTION ####

        ### GENERAL SERVER FUNCTIONS
        elif component_type =='server-helper-functions':
            function_type = data_options[0]
            if function_type == 'get_primary_geneID':
                geneID = get_primary_geneID_from_symbol_group(data_filter)
                #Send gene info data
                send_data(ws,'server-helper-functions','get_primary_geneID',json.dumps(gene_info))


##GENERAL SERVER FUNCTIONS
@app.route('/get_primary_gene_from_symbol_group', methods=['GET'])
def primary_gene_from_symbol_group():
    symbol_group = request.args.get('symbol_group')

    if not symbol_group:
        return jsonify({"error": "Missing 'symbol_group' parameter"}), 400

    primary_gene = get_primary_geneID_from_symbol_group(symbol_group)
    return jsonify({"primary_gene": primary_gene})

'''
    Helper function to get the primary GeneID from a symbol group
    @param symbol (string) => A symbol string
    @return geneID => Primary gene for that symbol group
'''
def get_primary_geneID_from_symbol_group(symbol):
    geneID = int( gene2pubmed_db[gene2pubmed_db['Symbol_group']==symbol]['Primary_gene'].value_counts().index[0])
    return geneID

#Function to generate summary_component_data
def generate_summary_component_data(db):
    """
        db = database to generate summary from
    """
    #Get Unique Publication IDs to filter gene2pubmed_papers_db
    publication_IDs = db['PubMed_ID'].drop_duplicates().values.tolist()
    paper_db = gene2pubmed_papers_db[gene2pubmed_papers_db['PubMed_ID'].isin(publication_IDs)]
    #Filter Research and Review Papers
    research_db = paper_db[paper_db['Publication_types'].str.contains('Research')]
    review_db = paper_db[paper_db['Publication_types'].str.contains('Review')]
    #Number of papers
    n_papers = len(paper_db)
    #Calculate value counts for histogram all paper types
    histogram_all, histogram_all_cum = calculate_histogram(paper_db['Pubdate'])
    #Number of Research papers
    n_research_papers = len(research_db)
    #Calculate value counts for histogram research papers 
    histogram_research, histogram_research_cum = calculate_histogram(research_db['Pubdate'])
    #Number of Review papers
    n_review_papers = len(review_db)
    #Calculate  value counts for histogram reviews papers
    histogram_reviews, histogram_reviews_cum = calculate_histogram(review_db['Pubdate'])

    #Rate Section
    histogram_all = calculate_acceleration_rate(histogram_all)
    histogram_research = calculate_acceleration_rate(histogram_research)
    histogram_reviews = calculate_acceleration_rate(histogram_reviews)
    return histogram_all, histogram_all_cum, n_papers, histogram_research, histogram_research_cum, n_research_papers, histogram_reviews, histogram_reviews_cum, n_review_papers



"""
    Helper function to calculate the acceleration rate
    @param hist (df) => list of values to calculate acceleration rate
    note: ID is the year for histogram afteer running calculate_histogram
"""
def calculate_acceleration_rate(hist):
    length = len(hist)
    acc_list = []
    #count to track index
    count = 0
    for i, each_row in hist.iterrows():
        if count != 0 and count != length:
            acc_rate = acceleration_rate(
                    each_row['counts'],
                    hist.iloc[count-1]['counts'],
                    each_row['ID'],
                    hist.iloc[count-1]['ID'])
            acc_list.append(acc_rate)
        else:
            acc_list.append(0)
        count = count + 1
    hist['acc']=acc_list
    return hist

"""
    Helper function to calculate the acceleration rate
    @param vel_end (int) => velocity at the end
    @param vel_start (int) => velocity at the start
    @param time_end (int) => year end
    @param time_start (int) => year start
"""
def acceleration_rate(vel_end,vel_start,time_end,time_start):
    vel_end = int(vel_end)
    vel_start = int(vel_start)
    time_end = int(time_end)
    time_start = int(time_start)
    try:
        acc = (vel_end - vel_start)/(time_end-time_start)
    except ZeroDivisionError:
        acc = 0
    return acc


"""
Helper function to transform geneID from NCBI to accession numbers.
Used now for adding the number to gene card. But this can allow us to connect to all the info in uniprot and possibly other dbs
"""
def ncbi_geneID_to_accession(geneID):
    try:    
        # Fetch the RefSeq record from the NCBI Nucleotide database for the gene
        handle =  Entrez.efetch(db="nuccore", id=geneID, rettype="fasta_cds_na", retmode="table")
        gene_record = handle.read()
        handle.close()
    except HTTPError:
        print('Error in response') 
        gene_record = ''
    # Extract the accession number from the gene record
    accession = None
    #Get the ascension number
    match = re.search(r"UniProtKB/Swiss-Prot:(\S+)", gene_record)
    if match:
        accession = match.group(1).rstrip("]")
        print(f"Accession number: {accession}")
    else:
        print("Accession number not found in response.")
    return accession

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
def search(query,db='gene',sort='relevance',retmax='20'):
    Entrez.email = 'agonzamart@gmail.com'
    handle = Entrez.esearch(db=db,
                            sort=sort,
                            retmax=retmax,
                            retmode='xml',
                            term=query)
    results = Entrez.read(handle)
    return results

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
def fetch_details(id_list,db='gene'):
    ids = ','.join(id_list)
    Entrez.email = 'agonzamart@gmail.com'
    handle = Entrez.efetch(db=db,
                           retmode='xml',
                           id=ids)
    results = Entrez.read(handle)
    # results2 = Entrez.parse(handle)
    return results

def create_gene_info_object(gene_data,GeneID):
    """
        Function to create gene_info_object to send to client
        with the parameters needed for gene_info page
        @params{
            gene_data: Data for a single Gene we want to parse.
    """
    #For the filters we need GeneID as an int
    GeneID = int(GeneID)
    #Create processing for what I want
    gene_info = {
        'Official_Symbol':'',
        'Official_Name':'',
        'Accession_Number':'',
        'Synonyms':'',
        'Type':'',
        'Organism':'',
        'Lineage':'',
        'Summary':'',
        'First_Published':'',
        'Last Published':'',
        'Total_Publications':'',
        'Peak_Publication_Year':'',
        'N_Publications_Peak_Year':'',
        'Avg_Publications_Year':'',
        'Total_Citations':'',
        'Peak_Citations_Year':'',
        'N_Publications_Peak_Year':'',
        'Avg_Citations_Year':'',
        'Peak_Citation_Year':'',
        'N_Citations_Peak_Year':'',
        'Total_Citations':'',
        'First_Citation': '',
        'Last_Citation':'',
        'Symbols_In_Group':''
    }
    #VALUES GATHERED FROM NCBI API
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
     #Accession Number
    try:
        gene_info['Accession_Number'] = ncbi_geneID_to_accession(GeneID)
    except KeyError:
        gene_info['Accession_Number'] = ''
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
    #### VALUES GATHERED FROM NCBI API END ####
    #Total Publications (Need to delete duplications that might occurr due to using primary gene)
    gene_info['Total_Publications'] = str(gene2pubmed_db[gene2pubmed_db['Primary_gene']==GeneID][['PubMed_ID','Pubdate']].drop_duplicates(subset='PubMed_ID').value_counts().sum())
    #First Published
    gene_info['First_Published'] = str(gene2pubmed_db[gene2pubmed_db['Primary_gene']==GeneID]['Pubdate'].min())
    #Last Published
    gene_info['Last_Published'] = str(gene2pubmed_db[gene2pubmed_db['Primary_gene']==GeneID]['Pubdate'].max())
    #Peak_Publication_Year, N_Publications_Peak_Year
    gene_info['Peak_Publication_Year'], gene_info['N_Publications_Peak_Year'] = gene_publication_counts(GeneID)
    #Avg_Publications_Year
    gene_info['Avg_Publications_Year'] = str(int(gene2pubmed_db[gene2pubmed_db['Primary_gene']==GeneID]['Pubdate'].value_counts().mean()))
    #Peak_Citation_Year, N_Citations_Peak_Year, Total_Citations
    gene_info['Peak_Citation_Year'], gene_info['N_Citations_Peak_Year'], gene_info['Total_Citations'] ,gene_info['Avg_Citations_Year'], gene_info['First_Citation'], gene_info['Last_Citation'] = gene_citation_counts(GeneID)
    gene_symbol = get_symbol_from_geneID(GeneID)
    gene_info['Symbols_In_Group'] = get_symbols_in_group(gene_symbol)
    return gene_info

"""
    Helper function that returns a symbol from a GeneID
    @param: GeneID
    @return: symbol name
"""
def get_symbol_from_geneID(GeneID):
    try:
        symbol = gene_info_db[gene_info_db['GeneID']==GeneID]['Symbol'].values.tolist()[0] 
    except IndexError:
        symbol = None
        print(str(GeneID) + ' id_not_found')
    return symbol

"""
    Helper function that returns all the unique symbols inside a symbol group
    @param: symbol_name: symnol name as string
    @return: symbol_in_group: list of strings
"""
def get_symbols_in_group(symbol_name):
    symbols_in_group = []
    genes_in_group =  symbol_group[symbol_group['Symbol']==symbol_name]['idList'].values.tolist()
    #Check if we have other genes in this group
    if len(genes_in_group)>0:
        for each_gene_in_group in symbol_group[symbol_group['Symbol']==symbol_name]['idList'].values.tolist()[0]:
            symbols_in_group.append(get_symbol_from_geneID(int(each_gene_in_group)))
    #If not is empty
    else:
        symbols_in_group = []
    symbols_in_group = sorted(set(symbol for symbol in symbols_in_group if symbol is not None))
    return symbols_in_group


"""
    Helper function to calculate the citation counts for a gene
    @param geneID: The geneID to search on
    @return top_year (year of max citations), citations_top_year (total citations year of max citations), total_citations (total amount of citations)
"""
def gene_citation_counts(geneID):
    #Filtering the database this is already in the code
    pubmed_ids = gene2pubmed_db[gene2pubmed_db['Primary_gene']==geneID]['PubMed_ID'].values.tolist()
    #Get the papers
    paper_db = gene2pubmed_papers_db[gene2pubmed_papers_db['PubMed_ID'].isin(pubmed_ids)]
    #This gets all the citations by year (We can use this to add a barchart of the citations every year in the gene card
    #Or somewhere else in the gene page
    gene_citations_by_year = paper_db.groupby('Pubdate')['Citations'].sum().reset_index().sort_values('Citations',ascending=False)
    total_citations = gene_citations_by_year['Citations'].sum()
    avg_citations = int(gene_citations_by_year['Citations'].mean())
    #We ordered by citations so top citation is first row
    gene_citations_top_year = gene_citations_by_year.iloc[0]
    top_year = gene_citations_top_year['Pubdate']
    citations_top_year = gene_citations_top_year['Citations']
    first_citation =  gene_citations_by_year['Pubdate'].min()
    last_citation = gene_citations_by_year['Pubdate'].max()
    return str(top_year), str(citations_top_year), str(total_citations), str(avg_citations), str(first_citation), str(last_citation)

"""
    @param: GeneID
    @return: max_year, count_value
"""
def gene_publication_counts(GeneID):
    #We need unique publications
    #Since publications can have multiple genes associated with them and we want to count each publication only once
    gene_publications =  gene2pubmed_db[gene2pubmed_db['Primary_gene']==GeneID].drop_duplicates(subset='PubMed_ID')
    gene_pub_counts = gene_publications['Pubdate'].value_counts()
    max_year = gene_pub_counts.idxmax()
    n_pub = gene_pub_counts[max_year]
    return str(max_year), str(n_pub)

def getCitationCount(publicationIDArray):
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
    hist = hist.to_frame(name="counts")
    #We need to add missing years
    #For that we create an additional DataFrame with all years and merge them
    hist_all_years = pd.DataFrame()
    if len(hist!=0):
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
    #If no data in histogram for the range we return an empty dataframe
    else:
        hist_all_years['ID'] = []
        hist_all_years['counts'] = []
        return hist_all_years,hist_all_years

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
            - index_col: Column used as reference for the counts
        Output: data_frame with column counts
    """
    counts_dataFrame = data[index_col].value_counts().to_frame(name='counts')
    #Symbols
    if index_col == 'Symbol_group':
        index_col='Symbol'
        #Reset index twice. First we extract the symbols and then we extract the ranking
        counts_dataFrame = counts_dataFrame.reset_index()
        counts_dataFrame = counts_dataFrame.reset_index()
        #Rename columns
        counts_dataFrame.rename(columns={"index": 'Ranking', "Symbol_group":index_col},inplace=True)
    #authors
    else:
        #flatten multiindex
        counts_dataFrame.index = [' '.join(a) for a in counts_dataFrame.index.to_flat_index()]
        index_col='Author'
        #Reset index twice. First we extract the symbols and then we extract the ranking
        counts_dataFrame = counts_dataFrame.reset_index()
        counts_dataFrame = counts_dataFrame.reset_index()
        #Rename columns
        counts_dataFrame.rename(columns={"level_0": 'Ranking', "index":index_col},inplace=True)
    column_order = ['Ranking',index_col,'counts']
    #We add 1 to start at 1 instead of at 0.
    counts_dataFrame['Ranking'] = counts_dataFrame['Ranking'] + 1
    #Reorder columnss
    counts_dataFrame = counts_dataFrame[column_order]
    return counts_dataFrame[:num_results]

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
    
