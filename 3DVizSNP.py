#!/usr/local/bin/python3

'''
3DVizSNP.py
Reads in a VCF file, submits variants to VEP server, generates an iCn3D link that 
shows the variants in the sequence track and highlights the mutant in the 3D viewer

Original script: Shashi Ratnayake, CGBB, CBIIT, NCI
Modifications by: Michael Sierk, NCI 
                  Manoj M Wagle, Université Grenoble Alpes; University of Sydney

TODO (5/16/23):
    - make gene selection more efficient by creating a bed file & intersecting with vcf
    - allow upload of a file with a list of gene ids
    - option to include SIFT/Polyphen score cutoff
    - load SIFT/PolyPhen scores into iCn3D
        - requires BED file, has to be loaded manually into iCn3D?
    - add other options for prediction, such as open cravat
'''

from argparse import ArgumentParser, HelpFormatter
import textwrap
import sys
import os.path
import re
from pysam import VariantFile
import requests # verify=False set in requests to avoid errors when on VPN
requests.packages.urllib3.disable_warnings() # removes InsecureRequestWarning due to verify=False being set
import json
import time
from datetime import datetime
from urllib.parse import quote 
#from urllib.request import urlopen 
import pandas as pd

    
def get_protein_id(gene, restapi):
    """ get Uniprot ID using Ensembl gene ID
    https://www.biostars.org/p/9529129/#9529154
    """

    SwissProt_ID = None

    # use REST API (default; slower, sometimes is down)
    if restapi:
        URL = 'https://rest.uniprot.org/idmapping'

        params = {
            'from': 'Ensembl',
            #'to': 'UniProtKB',
            'to': 'UniProtKB-Swiss-Prot',
            'ids': gene
        }

        response = requests.post(f'{URL}/run', params, verify=verify)

        job_id = response.json()['jobId']
        #print('job id:', job_id)
        job_status = requests.get(f'{URL}/status/{job_id}', verify=verify)
        d = job_status.json()

        # Make three attemps to get the results
        for i in range(3):
            #print(d.get('jobStatus'))
            if d.get("jobStatus") == 'FINISHED' or d.get('results'):
                job_results = requests.get(f'{URL}/results/{job_id}', verify=verify)
                results = job_results.json()
                #print(json.dumps(results, indent=2))
                for obj in results['results']:
                    SwissProt_ID = obj["to"]
                break
            time.sleep(1)
    else: 
        # can set up a local file to retrieve Ensembl->SwissProt mapping
        # use mapping file: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
        # Use following shell commands (can put into a shell script and run it):
        #  grep Ensembl HUMAN_9606_idmapping.dat| grep -v "_" > tmp.dat
        #  cut -f3 tmp.dat | cut -d"." -f1 > ens
        #  cut -f1 tmp.dat > sp
        #  paste ens sp > HUMAN_9606_idmapping_Ensembl.dat
        #  rm tmp.dat ens sp

        # creates a two column file:
        #  % head HUMAN_9606_idmapping_Ensembl.dat
        #  ENSG00000126653	A0A024QZ33
        #  ENSG00000249915	A0A024QZ42
        #  ENSG00000170312	A0A024QZP7
        # TODO: Need to deal with multiple UniprotIDs for 1 Ensembl ID.
        with open('HUMAN_9606_idmapping_Ensembl.dat') as f:
           ens_to_swissprot = dict([line.split() for line in f])
        SwissProt_ID = ens_to_swissprot[gene]

    return SwissProt_ID

def vep_output(variants, colnames):
    """ Run VEP with the identified variants and capture sift and polyphen scores"""
    
    global results

    species = args.s
    server = "https://rest.ensembl.org"
    ext = "/vep/" + species + "/region?uniprot=1"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    # combine all variants in a string to submit to VEP
    vcf_lines = "{\"variants\" : [" + variants + "]}"
    #print(vcf_lines)

    r = requests.post(server + ext, headers=headers, data=vcf_lines, verify=verify) 
    if not r.ok:
        r.raise_for_status()
        sys.exit()

    decoded = r.json()
    #with open('data.json', 'w') as f:
    #    json.dump(decoded, f)
    #print(json.dumps(decoded, indent=2))
    '''
    {
    "input": "1 3625748 rs200029021 T G .",
    "transcript_consequences": [
      {
        "transcript_id": "ENST00000344579",
        "consequence_terms": [
          "intron_variant"
        ],
      },
      {
        "codons": "gTg/gGg",
        "biotype": "protein_coding",
        "impact": "MODERATE",
        "polyphen_prediction": "probably_damaging",
        "cds_end": 329,
        "consequence_terms": [
          "missense_variant"
        ],
        "strand": 1,
        "gene_symbol_source": "HGNC",
        "variant_allele": "G",
        "cds_start": 329,
        "hgnc_id": "HGNC:27007",
        "cdna_start": 387,
        "protein_start": 110,
        "amino_acids": "V/G",
        "sift_score": 0,
        "protein_end": 110,
        "gene_id": "ENSG00000158109",
        "cdna_end": 387,
        "polyphen_score": 0.982,
        "uniparc": [
          "UPI000014067B"
        ],
        "swissprot": [
          "Q5T0D9.129"
        ],
        "transcript_id": "ENST00000378344",
        "gene_symbol": "TPRG1L",
        "sift_prediction": "deleterious"
      },
    '''

    vep_results = pd.DataFrame(columns=colnames)
    
    for i in decoded:

        var = i.get("input")

        for key in i:
            if key == "transcript_consequences":
                for j in i[key]:
                    # Note: retrieve all predictions, not just deleterious
                    if "sift_prediction" in j or "polyphen_prediction" in j:
                        gene_id = j.get("gene_id")
                        gene_sym = j.get("gene_symbol")

                        if "swissprot" in j:
                            sp = j.get("swissprot")[0].split('.')[0]

                            # deal with trembl-only entries later
                            #elif "trembl" in j:
                            #    sp = j.get("trembl")[0].split('.')[0]
                            #else:
                            #    sp = 'NA'

                            # get the amino acid mutation
                            aa = j.get("amino_acids").split("/")[0]  # A/C
                            aa = aa.strip('\s+')
                            aanum = str(j.get("protein_start"))
                            alt_aa = j.get("amino_acids").split("/")[1]
                            alt_aa = alt_aa.strip()
                            mutaa = aa + aanum + alt_aa  # e.g. T355C

                            row = pd.Series({'variant': var,
                                            'EnsID': gene_id,
                                            'Symbol': gene_sym,
                                            'SPID': sp,
                                            'PDBID': '',
                                            'mutaa': mutaa,
                                            'SIpred': j.get("sift_prediction", "NA"),
                                            'SIscore': j.get("sift_score", 99.9),
                                            'PPpred': j.get("polyphen_prediction", "NA"),
                                            'PPscore': j.get("polyphen_score", -1.0)})

                            #print("var:", var, "mutaa:", mutaa)
                            vep_results = pd.concat([vep_results, row.to_frame().T])

                            break # just take the first hit (avoid isoforms)

    vep_results.set_index('variant', inplace=True)

    return(vep_results)

def get_pdb_num(pdbid, spid, aanum):
    '''check for offset between Uniprot -> PDB residue number mapping

    "P01116": {
        "identifier": "RASK_HUMAN", 
        "mappings": [
          {
            "entity_id": 2, 
            "end": {
              "author_residue_number": null, 
              "author_insertion_code": "", 
              "residue_number": 190
            }, 
            "chain_id": "A", 
            "start": {
              "author_residue_number": 1, 
              "author_insertion_code": "", 
              "residue_number": 2
            }, 
            "unp_end": 189, 
            "unp_start": 1, 
            "struct_asym_id": "B"
          }
        ], 
        "name": "RASK_HUMAN"
      }
    '''
    if args.test:
        print("\n\t\tChecking numbering offset between UniProt", spid, "and PDBID", pdbid, "...")
        #print("\t\tpdb:", pdbid, "spid:", spid, "aanum:", aanum)

    pdb = pdbid.split("_")[0]
    chainid = pdbid.split("_")[1]
    pdb_url = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/" + pdb
    r = requests.get(pdb_url, verify=verify)
    #print(r.json())
    if r.status_code == 200:
        uniprot_map = r.json()
    else:
        uniprot_map = "None"
        return("None", "None")

    '''
                # from the RASK example above:
                #
                # unp_start = uniprot start #
                # author_residue_number = pdb # start
                # residue_number = position in pdb chain
                # e.g. KRAS 1VVB:
                # 0 1 2 3 4 5 6 7 8 9 10 11 12 13
                # G M T E Y K L V V V G  A  G*  G (* is mutated)
                #
                # mutaa = G12V, scap = 1VVB_A_12_V, sift = 13 V
                #
                # "unp_start": 1, M is the 1st residue of Uniprot P01116
                # "author_residue_number": 1, M is residue #1 in 1VVB (PDB start residue #)
                # "residue_number": 2, M is the 2nd position in the chain (starts with G0)
                #   (position in chain)
                #  -> what we need for iCn3D SIFT/PolyPhen track
                # so pdbnum = aanum - difference between uniprot and pdb numbering
                #           = aanum - (unp_start - author_residue_number)
                #    respos (for SIFT/PolyPhen in iCn3D) = 
                #           aanum - (uniprot_start - residue_number)
    '''        
    for id in uniprot_map:
        if id.upper() == pdb:
            mappings = uniprot_map[id]["UniProt"][spid]["mappings"]
            #mappings: [{'entity_id': 2, 'chain_id': 'E', 'struct_asym_id': 'C', 'start': {'author_residue_number': 6, 'author_insertion_code': '', 'residue_number': 1}, 'end': {'author_residue_number': 26, 'author_insertion_code': '', 'residue_number': 21}, 'unp_start': 111, 'unp_end': 131}]
            for i in mappings:
                if i['chain_id'] == chainid: 
                    uniprot_start_num = i['unp_start']
                    pdb_start_num = uniprot_start_num

                    if i['start']['author_residue_number'] != None:
                        pdb_start_num = i['start']['author_residue_number']
                    else:
                        pdb_start_num = i['start']['residue_number']
                    residue_position = i['start']['residue_number']
                    #print("uniprot_start_num:", uniprot_start_num, "pdb_start_num:", pdb_start_num, "residue_position:", residue_position)
                    
                    diff = int(uniprot_start_num) - int(pdb_start_num)

                    pdbnum = int(aanum) - diff 

                    # residue position in sequence needed for iCn3D tracks
                    respos = int(aanum) - (uniprot_start_num - residue_position)
                    respos = int(respos)

                    return(pdbnum, respos)
    
    return("None", "None") # if didn't match pdbid or chain_id 


def check_mutant(pdbid, chainid, pdbnum):

    if args.test:
        print("\t\tChecking if engineered mutant...")
    
    ''' Determine if amino acid is an engineered mutant in PDB structure

        https://www.ebi.ac.uk/pdbe/api/pdb/entry/mutated_AA_or_NA/1bgj
        {
            "1bgj": [
            {
                "entity_id": 1, 
                "residue_number": 116, 
                "author_residue_number": 116, 
                "chain_id": "A", 
                "author_insertion_code": "", 
                "mutation_details": {
                "to": "S", 
                "from": "C", 
                "type": "Engineered mutation"
            }, 
            "chem_comp_id": "SER", 
            "struct_asym_id": "A"
            }, 
            ...
            ]
        }
        '''

    mut = False

    url = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/mutated_AA_or_NA/" + pdbid 
    r = requests.get(url, verify=verify)
    #print(r.json())
    if r.status_code == 200:
        mut_list = r.json()
    else:
        mut_list = "None"
    
    if mut_list != "None":
        for m in mut_list[pdbid]:
            if m["chain_id"] == chainid and str(m["author_residue_number"]) == pdbnum:
                if m["mutation_details"]["type"] == "Engineered mutation":
                    mut = True

    return(mut)

def check_observed(pdbid, chainid, pdbnum):

    if args.test:
        print("\t\tChecking if residue is observed in structure...")
        #print("\t\t", pdbid, chainid, pdbnum)

    obs = True

    url = "https://www.ebi.ac.uk/pdbe/api/pdb/entry/residue_listing/" + pdbid + "/chain/" + chainid 
    r = requests.get(url, verify=verify)
    #print(r.json())
    # {'3n56': {'molecules': [{'entity_id': 2, 'chains': [{'struct_asym_id': 'C', 'residues': 
    # [{'residue_number': 1, 'residue_name': 'SER', 'author_residue_number': 1, 'author_insertion_code': '', 'observed_ratio': 1.0}
    if r.status_code == 200:
        res_dict = r.json()
    else:
        res_dict = "None"

    # convert author_residue_number and observed_ratio to dictionary
    res = [d.get('author_residue_number', None) for d in res_dict[pdbid]['molecules'][0]['chains'][0]['residues']]
    obsval = [d.get('observed_ratio', None) for d in res_dict[pdbid]['molecules'][0]['chains'][0]['residues']]
    resnum = dict(zip(res, obsval))

    if int(pdbnum) in resnum:
        if float(resnum[int(pdbnum)]) == 0.0:
            obs = False
    else:
        if args.test:
            print("\t\t\t", pdbnum, "not in list of residues, assuming numbering mismatch, setting obs to False.")
        
        obs = False

    return(obs)

def get_pdb_id(results):

    if args.test:
        print("\n\tget_pdb_id...")

    '''
    Retrieve the list of PDB IDs for each Uniprot ID
    For each amino acid position:
        1. identify if there is an x-ray or EM structure
        2. if so, get the PDB ID for the highest resolution structure
        3. if not, get the PDB ID of any NMR structures available
        4. fix offset between UniProt and PDB numbering
        5. check that residue appears in the experimental structure
        6. check that the residue in question is not an engineered mutant 
    '''

    spid_list = list(set(results.SPID)) # unique set of SwissProt ids

    for spid in spid_list:
        #print("spid:", spid)
        url = "https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/" + spid 
        r = requests.get(url, verify=verify)
        #print(r.json())
        if r.status_code == 200:
            pdbid_list = r.json()
        else:
            pdbid_list = "None"

        if (pdbid_list != 'None'):
            aa_list = results.loc[results['SPID']==spid, "mutaa"]

            for aa in aa_list:
                #print("aa:", aa)
                aanum = aa[1:-1] 

                for m in pdbid_list: # m is swissprot ID
                    maxres = 10

                    for j in pdbid_list[m]:
                        #print('pdbid:', j['pdb_id'],'chain:',j['chain_id']) #,'resolution:',j['resolution'],'start:', j["unp_start"], 'end:', j['unp_end'])
                        
                        # if n in range of pdb, use pdb id instead of uniprot id
                        if ((int(aanum) >= j['unp_start']) & (int(aanum) <= j['unp_end'])):
                            if j["resolution"] is None:
                                res = 0
                            else:
                                res = j["resolution"]

                            if (res == 0) & (maxres < 10):
                                continue # NMR, but already have xray

                            elif ((res == 0) & (maxres == 10) | (0 < res < maxres)):
                                # either NMR and no xray, or xray with better resolution
                                #   note: (unnecessarily) replaces existing PDBID if only have NMR structures

                                # account for offset between PDB and Uniprot numbering
                                pdbid = j["pdb_id"].upper() + "_" + j["chain_id"]
                                num = get_pdb_num(pdbid, spid, aanum)
                                if (num is None or "None" in num):
                                    continue
                                pdbnum = num[0]
                                respos = num[1]

                                #print("\tpdbid:", pdbid, "spid:", spid, "aanum:", aanum, "pdbnum:", pdbnum, "residue position:", respos)

                                # check if residue is observed in the structure                                
                                if not check_observed(j["pdb_id"], j["chain_id"], pdbnum):
                                    if args.test:
                                        print("\t\t  -> Residue", aa, "doesn't appear in structure:", j['pdb_id'], j['chain_id'], "(pdbnum:", pdbnum, ")")
                                    
                                    continue

                                # check if residue is an engineered mutant in PDB structure
                                if check_mutant(j["pdb_id"], j["chain_id"], str(pdbnum)):
                                    if args.test:
                                        print("PDB engineered mutation at", pdbnum, j["pdb_id"])
                                    
                                    continue

                                results.loc[(results['SPID']==spid) & (results['mutaa']==aa),'PDBID'] = pdbid
                                results.loc[(results['SPID']==spid) & (results['mutaa']==aa),'pdbnum'] = str(pdbnum)
                                # only need this to make SIFT & PolyPhen tracks, delete before html/csv output:
                                results.loc[(results['SPID']==spid) & (results['mutaa']==aa),'respos'] = respos

                                if res > 0:
                                    maxres = res # if xray, reset maxres

def variant_string(mut_list):
    ''' 
    create a string that includes all the variants from one PDB ID or Uniprot ID for iCn3D
    '''
    variant_str = ""
    for mutaa in mut_list:
        s = re.split(r'(\d+)', mutaa) # need the number and new aa
        if variant_str == '':
            variant_str += s[1] + ' ' + s[2]
        else:
            variant_str += ',' + s[1] + ' ' + s[2]
 
    return variant_str

def get_iCn3D_path(results):
    '''
    generates the iCn3D path(s) based on the variants for a given gene
    TODO: need to modify to produce BED file to show SIFT/Polyphen score
    '''

    date = datetime.now()
    url_path='https://www.ncbi.nlm.nih.gov/Structure/icn3d/full.html?'

    for row in results.itertuples():

        if args.test:
            print("\tGetting URL for", row.SPID, end='')

        iCn3Durl = "No structure for " + row.SPID

        sift_str = ''
        poly_str = ''
        scap_str = ''
        url_query = ''
        url_command = ''

        mutaa = row.mutaa
        s = re.split(r'(\d+)', mutaa)  # need the number and new aa

        # check to see if we are using AlphaFold structure
        # print('row.PDBID:', row.PDBID)
        if row.PDBID == '': # or row.PDBID == 'NA':
            # check length -> sequences > 2700 aa are not in AF predictions
            length_query = "https://rest.uniprot.org/uniprotkb/" + row.SPID + "?format=tsv&fields=length"
            r = requests.get(length_query, verify=verify).text
            if int(r.split()[1]) > 2700:
                print(row.SPID, "length > 2700, no AlphaFold prediction")
                continue
            else:
                if args.test:
                    print(" (AlphaFold)")

                if row.SIpred != '':
                    sift_str += s[1] + ' ' + s[2]
                if row.PPpred != '':
                    poly_str += s[1] + ' ' + s[2]

                scap_str += row.SPID + '_A' + '_' + s[1] + "_" + s[2]  # e.g. P16860_A_113_Y
            
                url_query =  'afid=' + row.SPID + '&date=' + date.strftime("%Y%m%d") + '&v=3.12.7&command='
                
                url_command = 'view annotations; set annotation cdd; set view detailed view;  set thickness | stickrad 0.2'    
                url_command += '; add track | chainid ' + row.SPID + '_A' + ' | title SIFT_predict | text ' + sift_str
                url_command += '; add track | chainid ' + row.SPID + '_A' + ' | title PolyPhen_predict | text ' + poly_str
                url_command += '; scap interaction ' + scap_str

        else:
            if args.test:
                print(' (PDB', str(row.PDBID), ')')

            if row.SIpred != '':
                sift_str += str(row.respos) + ' ' + s[2]

            if row.PPpred != '':
                poly_str += str(row.respos) + ' ' + s[2]

            scap_str += str(row.PDBID) + "_" + str(row.pdbnum) + "_" + s[2] # e.g. 1HLZ_A_113_Y

            url_query =  'pdbid=' + row.PDBID.split("_")[0] + '&date=' + date.strftime("%Y%m%d") + '&v=3.12.7&command='

            url_command = 'view annotations; set annotation cdd; set view detailed view;  set thickness | stickrad 0.2'    
            url_command += '; add track | chainid ' + row.PDBID + ' | title SIFT_predict | text ' + sift_str
            url_command += '; add track | chainid ' + row.PDBID + ' | title PolyPhen_predict | text ' + poly_str
            url_command += '; scap interaction ' + scap_str

        url_command = quote(url_command) # encode the spaces for URL
        iCn3Durl = url_path + url_query + url_command

        results.loc[row.Index,'Link'] = iCn3Durl

def print_html(results):
    '''
    Write out the results dataframe to an HTML file with iCn3D links
    '''

    fout = args.v + "_output.html"
    f = open(fout, 'w')
  
    # html code 
    html = """<html>
    <head>
    <title>iCn3D links</title>
    <style>
    table, th, td {
        border: 1px solid black;
        border-collapse: collapse;
    }
    th, td {
        padding: 10px;
    }
    th {
        background-color: #D3D3D3;
    }
    tr {
        border-bottom: 1px solid #ddd;
    }
    </style>
    </head>
    <body>
    <h2>Click on a link to open the deleterious variants in iCn3D</h2>
     <table border='1'>"""
    
    html += "Input file: " + args.v + "<br>"
    if args.t:
        html += "SIFT & PolyPhen scores taken from TCGA VCF file.<br>"

    # output table
    html += "<tr><th>#</th><th>Variant</th><th>GeneID</th><th>Gene Symbol</th>\
             <th>UniprotID</th><th>mutaa</th>\
             <th>PDB ID</th><th>pdbnum</th>\
             <th>SIFT Call</th><th>SIFT Score</th>\
             <th>PolyPhen Call</th><th>PolyPhen Score</th>\
             <th>iCn3D link</th></tr>"
    
    n = 0
    for row in results.itertuples():

        # limit output to 1000 rows by default
        n += 1
        if n > args.n:
            print("\tHTML output limited to", str(args.n), "rows...")
            break

        iCn3Durl = "<a href=" + str(row.Link) + " target=\"_blank\">iCn3D link</a><br>"

        # put in link for EnsID, SPID, PDBID
        # EnsID: https://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000187634
        # SPID: https://www.uniprot.org/uniprotkb/Q5T0D9
        # PDBID: https://www.rcsb.org/structure/1HLZ
        ensid_url = "https://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=" + row.EnsID
        ensid_link = "<a href=" + ensid_url + " target=\"_blank\">" + row.EnsID + "</a><br>"

        if row.SPID != '':
            spid_url = "https://www.uniprot.org/uniprotkb/" + row.SPID
            spid_link = "<a href=" + spid_url + " target=\"_blank\">" + row.SPID + "</a><br>"
        else:
            spid_link = row.SPID

        if row.PDBID != '':
            pdbid_url = "https://www.rcsb.org/structure/" + row.PDBID.split("_")[0]
            pdbid_link = "<a href=" + pdbid_url + " target=\"_blank\">" + row.PDBID + "</a><br>"
        else:
            pdbid_link = row.PDBID 

        # highlight deleterious mutations with color text
        SIpred_td = "<td>" + row.SIpred + "</td>"
        if row.SIpred == 'deleterious':
            SIpred_td = "<td style=\"color:red\">" + row.SIpred + "</td>"
        elif row.SIpred == 'deleterious_low_confidence':
            SIpred_td = "<td style=\"color:orange\">" + row.SIpred + "</td>"

        PPpred_td = "<td>" + row.PPpred + "</td>"
        if row.PPpred == 'probably_damaging':
            PPpred_td = "<td style=\"color:red\">" + row.PPpred + "</td>"
        elif row.PPpred == 'possibly_damaging':
            PPpred_td = "<td style=\"color:orange\">" + row.PPpred + "</td>"

        html += "<tr><td>" + str(n) + "</td><td>" + str(row.Index) + "</td><td>" \
            + ensid_link + "</td><td>" + row.Symbol + "</td><td>" \
            + spid_link + "</td><td>" + row.mutaa + "</td><td>" \
            + pdbid_link + "</td><td>" + str(row.pdbnum) + "</td>" \
            + SIpred_td + "<td>" + str(row.SIscore) + "</td>" \
            + PPpred_td + "<td>" + str(row.PPscore) + "</td><td>" \
            + iCn3Durl + "</td></tr>\n"

    html += "</table>"
    html += "</body>"
    html += "</html>"

    f.write(html)
    f.close()
    
def print_csv(results):
    '''
    print the results dataframe in a .csv file
    '''
    
    fout = args.v + "_results.csv"
    
    # put in link for EnsID, SPID, PDBID
    # EnsID: https://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=ENSG00000187634
    # SPID: https://www.uniprot.org/uniprotkb/Q5T0D9
    # PDBID: https://www.rcsb.org/structure/1HLZ

    results_csv = results.copy()

    for row in results_csv.itertuples():

        ensid_link = '=HYPERLINK("https://useast.ensembl.org/Homo_sapiens/Gene/Summary?g=' + row.EnsID + '","' + row.EnsID + '")' 

        if row.SPID != '':
            spid_link = '=HYPERLINK("https://www.uniprot.org/uniprotkb/' + row.SPID + '","' + row.SPID + '")' 
        else:
            spid_link = row.SPID

        if row.PDBID != '':
            pdbid_link = '=HYPERLINK("https://www.rcsb.org/structure/' + row.PDBID.split("_")[0] + '","' + row.PDBID + '")' 
        else:
            pdbid_link = row.PDBID 
    
        iCn3Dlink = '=HYPERLINK("' + str(row.Link) + '","iCn3D Link")'

        results_csv.loc[row.Index,'EnsID'] = ensid_link
        results_csv.loc[row.Index,'SPID'] = spid_link
        results_csv.loc[row.Index,'PDBID'] = pdbid_link
        results_csv.loc[row.Index,'Link'] = iCn3Dlink

    results_csv.reset_index(inplace=True)

    csv_header = ['Variant', 'EnsemblID', 'GeneSymbol', 
                  'UniprotID', 'AminoAcidMutation', 
                  'PDBID', 'PDBnum',  
                  'SiftPrediction', 'SiftScore',
                  'PolyPhenPrediction', 'PolyPhenScore', 'iCn3Dlink']

    results_csv.to_csv(fout,  header=csv_header) # index=False, index_label='Variant',header=csv_header

def get_vcf():
    '''
    read in the VCF file & extract variants
    '''

    vcf = []
    num_lines = 0
    vcf_reader = VariantFile(args.v)

    '''
    convert gene symbols to ensembl gene ids:
 
    decoded = r.json()
    print(repr(decoded))
    "BRAF": {
        "version": 14,
        "seq_region_name": "7",
        "start": 140719327,
        "object_type": "Gene",
        "canonical_transcript": "ENST00000646891.2",
        "source": "ensembl_havana",
        "assembly_name": "GRCh38",
        "db_type": "core",
        "id": "ENSG00000157764",
        "strand": -1,
        "biotype": "protein_coding",
        "species": "homo_sapiens",
        "end": 140924929,
        "display_name": "BRAF",
        "description": "B-Raf proto-oncogene, serine/threonine kinase [Source:HGNC Symbol;Acc:HGNC:1097]",
        "logic_name": "ensembl_havana_gene_homo_sapiens"
    },    
    '''
    if args.g:
        #get gene chr and coordinates from ensembl
        genes = args.g.split(",")

        server = "https://rest.ensembl.org"
        ext = "/lookup/symbol/homo_sapiens"
        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
        #r = requests.post(server+ext, headers=headers, data='{ "symbols" : ["BRCA2", "BRAF" ] }')
        data = '{"symbols" :' + json.dumps(genes) + '}' # use json.dumps to get double quotes
        r = requests.post(server+ext, headers=headers, data=data, verify=verify)

        if not r.ok:
            r.raise_for_status()
            sys.exit()
        decoded = r.json()

        # for now just do one at a time 
        # (TODO create a bed file & intersect: https://www.biostars.org/p/46331/)
        for i in decoded:
            chr = "chr" + decoded[i]["seq_region_name"]
            start = decoded[i]["start"]
            stop = decoded[i]["end"]
            # print("i:", i)
            # print("display_name:", decoded[i]["display_name"])
            # print("id:", decoded[i]["id"])
            for record in vcf_reader.fetch(chr, int(start), int(stop+1)):
                num_lines += 1
                line = str(record)
                vcf_line = line.split()
                # only get variants
                if (vcf_line[3] != vcf_line[4]) and (vcf_line[4] != '.'):
                    vcf_line[0] = vcf_line[0].replace("chr","") # can't have chr in VEP REST API submission
                    vcf.append(vcf_line[0:6])

    else:
        for record in vcf_reader.fetch():
            num_lines += 1
            line = str(record)
            vcf_line = line.split()
            # only get variants
            if (vcf_line[3] != vcf_line[4]) and (vcf_line[4] != '.'):
                vcf_line[0] = vcf_line[0].replace("chr","") # can't have chr in VEP REST API submission
                vcf.append(vcf_line[0:6])

    return(vcf)

def get_vcf_tcga(colnames):
    ''' 
    Get Ensembl ID, SwissProt ID, SIFT, Polyphen predictions & scores directly from TCGA VCF file
    Note: not recommended, there are lots of discrepancies between TCGA values & VEP/Ensembl REST API values
    '''

    # read in the VCF file
    vcf = []
    vcf_reader = VariantFile(args.v)
    tcga_results = pd.DataFrame(columns=colnames)

    for record in vcf_reader.fetch():
        line = str(record)
        vcf_line = line.split()
        variant, mutaa, spid = '', '', ''
        c,d,e,sift_split,polyph_split = (), (), (), (), ()
        
        # only get variants
        if (vcf_line[3] != vcf_line[4]) and (vcf_line[4] != '.'):
            vcf_line[0] = vcf_line[0].replace("chr","") 
            variant = vcf_line[0:6]
            vcf.append(variant)

            c = record.info["CSQ"] # consequence: multiple INFO strings stored as a tuple
            #Gene: position 5
            #Protein: 15 110/272
            #Amino acids: 16 V/G
            #SwissProt: 31
            #SIFT: 36 deleterious(0)
            #PolyPhen: 37 probably_damaging(0.993)

            for d in c: 
                e = d.split("|")
                if e[1] == "missense_variant":
                    #print("gene id:", e[4], "swissprot:", e[31]) # lots of Trembl IDs, need Swissprot to get PDB IDs
                    if e[31] == '': # lot of these are missing in TCGA files
                        spid = get_protein_id(e[4], 1)
                        print("getting spid from REST API:", spid)
                    else:
                        spid = e[31]

                    mutaa = e[15].split("/")[0] + e[14].split("/")[0] + e[15].split("/")[1]

                    sift_split = ['','',0,'']
                    if e[35] != '':
                        sift_split = re.split(r'(\w+)\(([0-9]+\.?\d*)\)', e[35])

                    polyph_split = ['','',0,''] 
                    if e[36] != '':
                        polyph_split = re.split(r'(\w+)\(([0-9]+\.?\d*)\)', e[36])

                    row = pd.Series({'variant': variant,
                                    'EnsID': e[4],
                                    'Symbol': e[3],
                                    'SPID': spid,
                                    'PDBID': '',
                                    'mutaa': mutaa,
                                    'SIpred': sift_split[1],
                                    'SIscore': float(sift_split[2]),
                                    'PPpred': polyph_split[1],
                                    'PPscore': float(polyph_split[2])})
                    tcga_results = pd.concat([tcga_results, row.to_frame().T])

    tcga_results.set_index('variant', inplace=True)

    return(vcf, tcga_results)

def cli():

    # format the description
    class RawFormatter(HelpFormatter):
        def _fill_text(self, text, width, indent):
            return "\n".join([textwrap.fill(line, width) for line in textwrap.indent(textwrap.dedent(text), indent).splitlines()])

    desc = f'''
            SNP2iCn3D.py: Runs VEP on single-nucleotide variants extracted from a VCF file 
            and generates iCn3D links for predicted deleterious mutations.  
            If the VCF file is from TCGA, the -t flag extracts the VEP data from the VCF file instead of 
            submitting to the VEP server.

            Two modes: 
                1) provide a comma-separated list of gene symbols with the -g flag, 
                   only locations matching those genes will be extracted; 
                2) do not provide the -g flag, all variants 
                   in the VCF will be run through VEP. 
                    
            Output:    
                - an html file listing the variants, SIFT & PolyPhen scores, and iCn3D links
                - a .csv file that can be imported into Numbers or Google Sheets (Note that the iCn3D URLs are broken in Excel)
            '''
    def check_max(n):
        n = int(n)
        # restrict output lines to 10,000
        if (n > 10000):
            msg = "%r cannot be larger than 10,000, resetting value to 10,000." % n
            n = 10000
            raise ArgumentTypeError(msg)
        return n
    
    parser = ArgumentParser(description=desc, formatter_class=RawFormatter)

    parser.add_argument('-v', required=True, metavar='VCF', help="VCF file to extract the variants; must be compressed with bgzip \
        and the tabix .tbi index file must be present.")
    
    parser.add_argument('-n', type=check_max, default=1000, help="Number of output rows (default 1000, max 10,000)")

    parser.add_argument('-g', metavar='GENE', help="Select only variants from gene symbols of interest")
    
    parser.add_argument('-t', action='store_true', help="(deprecated) Extract SIFT & PolyPhen scores from VCF file from TCGA instead of submitting to VEP")
                        #action='store_true' means default is false 
    
    parser.add_argument('-s', type=str, default='human', help="species (default human) (use a common name from http://rest.ensembl.org/info/species.json)")
    # human (homo_sapiens)
    # mouse (mus_musculus)
    # rat (rattus_norvegicus)
    # dog (canis_lupus_familiaris)
    # chicken (gallus_gallus)
    
    parser.add_argument('--verifyfalse', action='store_false', 
                        help="If you are on a VPN and get error messages from the requests library, set the --verifyfalse flag")

    parser.add_argument('--test', action='store_true', help="Perform testing functions") # default is false

    return parser

def main():
    
    global args, verify
    args = cli().parse_args()
    verify = args.verifyfalse # default true

    # set up results dataframe
    colnames = ["variant", "EnsID", "Symbol", "SPID", "mutaa", "PDBID", "pdbnum", 
                "SIpred", "SIscore", "PPpred", "PPscore", "Link", "respos"]
    results = pd.DataFrame(columns=colnames)
    results.set_index('variant', inplace=True)

    def estimate_time(vcfdf):

        print("\n################# Time to Complete Job Estimate #################")

        l = len(vcfdf)
        total_batches = l/200
        time_per_batch = 10 # seconds
        mutations_per_batch = 12 # just an estimate 
        max_batches = int(args.n/mutations_per_batch) # how many batches to get to n (default 1000)

        if total_batches > max_batches:
            max_batches = max_batches
        else:
            max_batches = total_batches
        
        print("\tLength of vcf file:", l)
        print("\tEstimated maximum number of batches:", max_batches)
        
        vep_time = time_per_batch * max_batches
        pdb_time = 10 * min(args.n, l)/10 # ~10% of variants require PDB calls
        time = vep_time + pdb_time
        time = int(time/60 + 1)

        # Length of vcf file: 14
        # Estimated # of batches of 200 variants needed to get to 1000 mutations: 0.07
        # Estimated time at 10 seconds per batch if 12.5 mutations per batch: 17 minutes

        print("\tEstimated # of batches of 200 variants needed to get to", args.n, "mutations:", str(max_batches))
        print("\tEstimated time at", str(time_per_batch), "seconds per batch if", 
              str(mutations_per_batch), "mutations per batch:", str(time), "minutes")
        print("################# Time to Complete Job Estimate #################")

    if args.test and os.path.exists("vep_results.csv"):
        print("reading vep_results.csv...")
        # variant	EnsID	Symbol	SPID	mutaa	PDBID	pdbnum	SIpred	SIscore	PPpred	PPscore	Link
        # 1 943298 . C T .	ENSG00000187634	SAMD11	Q96NU1	P537L			tolerated	0.93	benign	0.009	
        results = pd.read_csv("vep_results.csv", index_col='variant') 
                                                 
        results['EnsID'] = results['EnsID'].astype("string")
        results['mutaa'] = results['mutaa'].astype("string")
        results['PDBID'] = results['PDBID'].astype("string")
        results['SIpred'] = results['SIpred'].astype("string")
        results['PPpred'] = results['PPpred'].astype("string")
        results['Link'] = results['Link'].astype("string")

    else:

        # Extract SIFT & PolyPhen scores from TCGA file
        if args.t:
            print("Getting VEP values from TCGA file...")
            vcf, results = get_vcf_tcga(colnames)  # returns vcf, results df

            if not vcf:
                sys.exit("No variants found in VCF file! Ending program...")

        # Get variants from VCF file, SIFT & PolyPhen from VEP REST API (default)
        else:
            # extract the variants from the vcf file
            vcf = get_vcf()
            if not vcf:
                sys.exit("No variants found in VCF file! Ending program...")
            else:
                estimate_time(vcf)

            # need to break into chunks of 200 variants - maximum POST size is 200
            def divide_variants_list(list, n):
                for i in range(0, len(list), n):
                    yield list[i:i + n]

            variant_list = list(divide_variants_list(vcf, 200))  # returns a list of lists

            # combine all variants in a string to submit to VEP
            print("\nSubmitting variants to VEP server (batch = 200 variants)...")
            n = 0
            for v in variant_list:
                n += 1
                print("\tprocessing batch", n, "...")
                variants = ''
                for c in v:
                    loc = " ".join(c)
                    variants += "\"" + loc + "\" ,"
                variants = variants[:-1]

                vep_result = []
                vep_result = vep_output(variants, colnames)

                results = pd.concat([results, vep_result])

                if len(results) >= args.n:
                    print("Stopping after", str(args.n), "missense mutations...")
                    break

                time.sleep(2)  # avoid taxing server
            print("Done")

            if args.test:
                results.to_csv("vep_results.csv")
        
        # end if args.t (TCGA)

    # end if args.test and os.path.exists("vep_results.csv"):

    if args.test:
        results['PDBID'].fillna('',inplace=True)  # need when testing and reading from .csv file
        results['SIpred'].fillna('',inplace=True) # need when testing and reading from .csv file
        results['PPpred'].fillna('',inplace=True) # need when testing and reading from .csv file

    # find if mutations are in PDB structures or not
    print("\nGetting PDB IDs...")
    print("\tChecking numbering offset, if residue is observed in structure, if residue is an engineered mutant...")

    get_pdb_id(results)
    print("Done") 

    print("\nGenerating iCn3D URLs...")
    get_iCn3D_path(results)

    # get rid of residue position (only need to make SIFT & PolyPhen tracks, don't want it in output table)
    results = results.drop(columns=['respos'])

    # generate an html page with results dataframe, iCn3D links
    print("Writing html file...")
    print_html(results)  

    # generate a .csv file with results dataframe, iCn3D links
    print("Writing .csv file...")
    print_csv(results)  

if __name__ == '__main__':   
    main()   
