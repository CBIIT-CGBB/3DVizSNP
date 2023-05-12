#!/usr/local/bin/python3

'''
3DVizSNP.py
Reads in a VCF file, submits variants to VEP server, generates an iCn3D link that 
shows the variants in the sequence track and highlights the first mutant in the 3D viewer

Original script: Shashi Ratnayake, CGBB, CBIIT, NCI
Modifications by: Michael Sierk, NCI 
                  Manoj M Wagle, Université Grenoble Alpes; University of Sydney

TODO (3/13/23):
    - make requests verify=False a cli flag
    - fix links so only 1 mutation per link (?)
    - add option to submit HUGO gene names instead of Ensembl IDs
    - option to include SIFT/Polyphen score cutoff
        - look at other SIFT/Polyphen classifications (e.g. Possibly Damaging)
    - load SIFT/PolyPhen scores into iCn3D
        - requires BED file, has to be loaded manually into iCn3D?
    - add other options for prediction, such as open cravat
'''
from argparse import ArgumentParser, HelpFormatter
import textwrap
from collections import defaultdict
import sys
import re
import re
from pysam import VariantFile
import requests # verify=False set in requests to avoid errors when on VPN
requests.packages.urllib3.disable_warnings() # removes InsecureRequestWarning due to verify=False being set
#import json
import time
from datetime import datetime
from urllib.parse import quote 
from urllib.request import urlopen 
import pandas as pd


    
def get_protein_id(gene, rest):
    """ get Uniprot ID using Ensembl gene ID
    https://www.biostars.org/p/9529129/#9529154
    """
    #print("getting Uniprot ID for", gene)
    SwissProt_ID = None

    # use REST API (default; slower, sometimes is down)
    if rest:
        URL = 'https://rest.uniprot.org/idmapping'

        params = {
            'from': 'Ensembl',
            #'to': 'UniProtKB',
            'to': 'UniProtKB-Swiss-Prot',
            'ids': gene
        }

        response = requests.post(f'{URL}/run', params, verify=False)
        #print(response)

        job_id = response.json()['jobId']
        #print('job id:', job_id)
        job_status = requests.get(f'{URL}/status/{job_id}', verify=False)
        d = job_status.json()

        # Make three attemps to get the results
        for i in range(3):
            #print(d.get('jobStatus'))
            if d.get("jobStatus") == 'FINISHED' or d.get('results'):
                job_results = requests.get(f'{URL}/results/{job_id}', verify=False)
                results = job_results.json()
                #print(json.dumps(results, indent=2))
                for obj in results['results']:
                    SwissProt_ID = obj["to"]
                break
            time.sleep(1)
    else: 
        # set up a local file to retrieve Ensembl->SwissProt mapping
        # can use mapping file: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz
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
        #ENSG00000131686	P23280
        #ENSG00000131686	Q8N4G4

    return SwissProt_ID

def vep_output(args, variants, colnames):
    """ Run VEP with the identified variants and capture sift and polyphen scores"""
    
    global results

    species = args.s
    print("species:", species)
    server = "https://rest.ensembl.org"
    ext = "/vep/" + species + "/region?uniprot=1"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}

    # combine all variants in a string to submit to VEP
    vcf_lines = "{\"variants\" : [" + variants + "]}"
    #print(vcf_lines)

    r = requests.post(server + ext, headers=headers, data=vcf_lines, verify=False) 
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
        #print("input:", var, "\n")
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
    #print(vep_results)

    return(vep_results)

def get_pdb_id(results):
    '''
    Retrieve the list of PDB IDs for each Uniprot ID
    For each amino acid position:
        1. identify if there is an x-ray or EM structure
        2. if so, get the PDB ID for the highest resolution structure
        3. if not, get the PDB ID of any NMR structures available
    '''
    spid_list = list(set(results.SPID)) # unique set of SwissProt ids

    for spid in spid_list:
        #print("spid:", spid)
        url = "https://www.ebi.ac.uk/pdbe/graph-api/mappings/best_structures/" + spid 
        r = requests.get(url, verify=False)
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
                for m in pdbid_list:
                    maxres = 10
                    for j in pdbid_list[m]:
                        #print(aanum)
                        #print('pdbid:', j['pdb_id'],'chain:',j['chain_id'],'resolution:',j['resolution'],'start:', j["unp_start"], 'end:', j['unp_end'])
                        # if n in range of pdb, use pdb id instead of uniprot id
                        if ((int(aanum) >= j['unp_start']) & (int(aanum) <= j['unp_end'])):
                            #print(type(j['resolution']),' resolution:',j['resolution'],"|",sep='',)
                            if j["resolution"] is None:
                                res = 0
                            else:
                                res = j["resolution"]
                            #print("res:", res, "maxres:", maxres)
                            if (res == 0) & (maxres < 10):
                                # NMR, but already have xray
                                break
                            elif ((res == 0) & (maxres == 10) | (0 < res < maxres)):
                                # either NMR, no xray, or xray with better resolution
                                # note: replaces existing PDBID if NMR only
                                pdbid = j["pdb_id"].upper() + "_" + j["chain_id"]
                                results.loc[(results['SPID']==spid) & (results['mutaa']==aa),'PDBID'] = pdbid
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
    
    Currently generates a separate URL for each PDB/AFID, but probably need
    to go back to a single URL for each mutation, as it is confusing when opening
    multiple mutations at a time in iCn3D
    '''
    
    date = datetime.now()
    url_path='https://www.ncbi.nlm.nih.gov/Structure/icn3d/full.html?'

    for row in results.itertuples():
        print("\nGetting URL for", row.SPID)
        iCn3Durl = "No structure for " + row.SPID

        sift_str = ''
        poly_str = ''
        scap_str = ''
        url_query = ''
        url_command = ''

        mutaa = row.mutaa
        s = re.split(r'(\d+)', mutaa)  # need the number and new aa

        if row.SIpred != '':
            sift_str += s[1] + ' ' + s[2]

        if row.PPpred != '':
            poly_str += s[1] + ' ' + s[2]

        # check to see if we are using AlphaFold structure
        if (row.PDBID == ''):
            # check length -> sequences > 2700 aa are not in AF predictions
            length_query = "https://rest.uniprot.org/uniprotkb/" + row.SPID + "?format=tsv&fields=length"
            r = requests.get(length_query, verify=False).text
            if int(r.split()[1]) > 2700:
                print(row.SPID)
                print("length > 2700, no AF prediction")
                continue
            else:
                print("Getting alphafold url...")
                scap_str += row.SPID + '_A' + '_' + s[1] + "_" + s[2]  # e.g. P16860_A_113_Y
            
                url_query =  'afid=' + row.SPID + '&date=' + date.strftime("%Y%m%d") + '&v=3.12.7&command='
                
                url_command = 'view annotations; set annotation cdd; set view detailed view;  set thickness | stickrad 0.2'    
                url_command += '; add track | chainid ' + row.SPID + '_A' + ' | title SIFT_predict | text ' + sift_str
                url_command += '; add track | chainid ' + row.SPID + '_A' + ' | title PolyPhen_predict | text ' + poly_str
                url_command += '; scap interaction ' + scap_str

        else:
            print('getting PDB url...')

            # have to check for offset between Uniprot -> PDB residue number mapping
            pdb = row.PDBID.split("_")[0]
            chainid = row.PDBID.split("_")[1]
            pdb_url = "https://www.ebi.ac.uk/pdbe/api/mappings/uniprot_segments/" + pdb
            r = requests.get(pdb_url, verify=False)
            if r.status_code == 200:
                uniprot_map = r.json()
            else:
                uniprot_map = "None"
        
            for id in uniprot_map:
                if id.upper() == pdb:
                    mappings = uniprot_map[id]["UniProt"][row.SPID]["mappings"]
                    if mappings[0]['chain_id'] == chainid:
                        pdb_start_num = mappings[0]['start']['author_residue_number']
                        uniprot_start_num = mappings[0]['unp_start']
        
            diff = 0
            if pdb_start_num != None:
                diff = uniprot_start_num - int(pdb_start_num)

            pdb_mutaa_num = int(s[1]) - diff # correct for uniprot->pdb offset 
            scap_str += row.PDBID + "_" + str(pdb_mutaa_num) + "_" + s[2] # e.g. 1HLZ_A_113_Y

            url_query =  'pdbid=' + row.PDBID.split("_")[0] + '&date=' + date.strftime("%Y%m%d") + '&v=3.12.7&command='

            url_command = 'view annotations; set annotation cdd; set view detailed view;  set thickness | stickrad 0.2'    
            url_command += '; add track | chainid ' + row.PDBID + ' | title SIFT_predict | text ' + sift_str
            url_command += '; add track | chainid ' + row.PDBID + ' | title PolyPhen_predict | text ' + poly_str
            url_command += '; scap interaction ' + scap_str

        url_command = quote(url_command) # encode the spaces for URL
        iCn3Durl = url_path + url_query + url_command

        results.loc[row.Index,'Link'] = iCn3Durl

def print_html(args, results):
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
    html += "<tr><th>#</th><th>Variant</th><th>GeneID</th><th>Symbol</th><th>UniprotID</th><th>PDB ID</th>\
             <th>mutaa</th><th>SIFT Call</th><th>SIFT Score</th><th>PolyPhen Call</th><th>PolyPhen Score</th>\
             <th>iCn3D link</th></tr>"
    
    n = 0
    for row in results.itertuples():

        # limit output to 1000 rows by default
        n += 1
        if n > args.n:
            print("HTML output limited to", str(args.n), "rows...")
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
            + ensid_link + "</td><td>" + row.Symbol + "</td><td>" + spid_link + "</td><td>" \
            + pdbid_link + "</td><td>" + row.mutaa + "</td>" \
            + SIpred_td + "<td>" + str(row.SIscore) + "</td>" \
            + PPpred_td + "<td>" + str(row.PPscore) + "</td><td>" \
            + iCn3Durl + "</td></tr>\n"

    html += "</table>"
    html += "</body>"
    html += "</html>"

    f.write(html)
    f.close()
    
def print_csv(args, results):
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
            pdbid_link = '=HYPERLINK("https://www.rcsb.org/structure/' + row.PDBID + '","' + row.PDBID + '")' 
        else:
            pdbid_link = row.PDBID 
    
        iCn3Dlink = '=HYPERLINK("' + row.Link + '","iCn3D Link")'

        results_csv.loc[row.Index,'EnsID'] = ensid_link
        results_csv.loc[row.Index,'SPID'] = spid_link
        results_csv.loc[row.Index,'PDBID'] = pdbid_link
        results_csv.loc[row.Index,'Link'] = iCn3Dlink

    results_csv.reset_index(inplace=True)

    csv_header = ['Variant', 'EnsemblID', 'Symbol', 'UniprotID', 'PDBID', 'AminoAcidMutation', 
                  'SiftPrediction', 'SiftScore',
                  'PolyPhenPrediction', 'PolyPhenScore', 'iCn3Dlink']

    results_csv.to_csv(fout,  header=csv_header) # index=False, index_label='Variant',header=csv_header

def get_vcf(args):
    '''
    read in the VCF file & extract variants
    '''
    vcf = []
    vcf_reader = VariantFile(args.v)
    n = 0
    for record in vcf_reader.fetch():
        n += 1
        line = str(record)
        vcf_line = line.split()
        # only get variants
        if (vcf_line[3] != vcf_line[4]) and (vcf_line[4] != '.'):
            vcf_line[0] = vcf_line[0].replace("chr","") # can't have chr in VEP REST API submission
            vcf.append(vcf_line[0:6])

            if len(vcf) == args.n:
                print("Stopping after", str(args.n), "variants out of", n, "vcf lines)")
                break

    return(vcf)

def get_vcf_tcga(args, colnames):
    ''' 
    Get Ensembl ID, SwissProt ID, SIFT, Polyphen predictions & scores directly from TCGA VCF file
    Note: not recommended, there are lots of discrepancies between TCGA values & VEP/Ensembl REST values
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

            if len(vcf) == args.n:
                print("Stopping after", str(args.n), "variants...")
                break

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
                1) provide a comma-separated list of specific Ensembl genes with the -g flag, 
                   only locations matching those genes will be extracted; 
                2) do not provide the -g flag, all variants 
                   in the VCF will be run through VEP. 
                    
            Output:    
                - an html file listing the variants, SIFT & PolyPhen scores, and iCn3D links
                - a .csv file that can be imported into Numbers or Google Sheets (URLs are broken in Excel)
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
    parser.add_argument('-g', metavar='GENE', help="Select only variants from Ensembl Gene IDs of interest")
    
    parser.add_argument('-v', required=True, metavar='VCF', help="VCF file to extract the variants; must be compressed with bgzip \
        and the tabix .tbi index file must be present.")
    
    parser.add_argument('-t', action='store_true', help="Extract SIFT & PolyPhen scores from VCF file from TCGA instead of submitting to VEP")
                        #action='store_true' means default is false 
    
    parser.add_argument('-s', type=str, default='human', help="species (default human) (use a common name from http://rest.ensembl.org/info/species.json)")
    
    parser.add_argument('-n', type=check_max, default=1000, help="Number of output rows (default 1000, max 10,000)")

    return parser

def main(args):
    '''
    Output:
        Input VCF file: <vcf_file> # could name <vcf_file>_out.csv?
        Coord EnsGene SPid PDBID mutaa SIFT Polyphen iCn3Dlink
          -- mutaa, iCn3dlink grouped by pdbid, spid
    
     -read in the whole VCF file, use client.get_gene_ids to fill gene_ids dictionary: 
        gene_ids[loc]["ens_id"] and gene_ids[loc]["sp_id"]
        -if args.g, select subset of genes
     -go through genes, get VEP results 
        -if args.t, get EnsID, SwissProtID, SIFT, Polyphen results from TCGA VCF file (SwissProt ID may not be current -> replace with get_protein_id)
     -generate iCn3D link for each gene
    
     Data structures
      vcf: list of single nucleotide variants pulled from VCF file "19 15256965 . T G . . ."
      gene_ids: gene_ids[loc]["ens_id"] and gene_ids[loc]["sp_id"]
        -> requires vcf
        -> limited to genes in args.g if present
      gene_id_list: list of Ensembl gene ids
        -> extracted from gene_ids
      variants (subset of vcf): string of variants for submitting to vep "19 15256965 . T G . . ."
        -> requires vcf, gene_ids_select
      sift: dictionary sift["gene_id"][mutaa] = {"sift_prediction": j.get("sift_prediction"), "sift_score": j["sift_score"]}
        -> requires variants
      polyphen: dictonary polyphen["gene_id"][mutaa] = {"polyphen_prediction": j.get("polyphen_prediction"), "polyphen_score": j["polyphen_score"]}
        -> requires variants
      gene_to_pid: dictionary gene_to_pid['ens_id'] = spid
        -> requires gene_ids
     
     Pandas data frame including all of above:
      variant EnsID SPID mutaa PDBID SIpred SIscore PPpred PPscore 

      url_list: dictionary url_list[gene] 
        -> requires sift, polyphen, gene, gene_to_pid[gene]
        -> 1..n urls per gene
        iCn3D link:
            - each gene has a list of URLs (1...n)
            - need afid/pdbid for each aa position
                struct_id[n] = pid or pdbid
            - need SIFT, PolyPhen strings for each afid/pdbid 
                variant_string
            - need scap strings for each afid/pdbid (different from SIFT/PolyPhen strings)

     HTML:
     EnsGene SPid PDBID SIFT PolyPhen iCn3Dlink
      - all mutations for a given PDB/SPid together in same iCn3D link
        variant_string2
     .csv:
     Coord EnsGene SPid PDBID mutaa SIFT Polyphen iCn3Dlink
    '''

    colnames = ["variant", "EnsID", "Symbol", "SPID", "PDBID", "mutaa", 
                "SIpred", "SIscore", "PPpred", "PPscore", "Link"]
    results = pd.DataFrame(columns=colnames)
    results.set_index('variant', inplace=True)

    # Extract SIFT & PolyPhen scores from TCGA file
    if args.t:
        print("Getting VEP values from TCGA file...")
        vcf, results = get_vcf_tcga(args, colnames) # returns vcf, gene_ids dict, fills sift & polyphen
        
        # find if mutations are in PDB structures or not
        print("\nGetting PDB IDs...", end='')
        get_pdb_id(results)
        print("Done") 

    # Get variants from VCF file, SIFT & PolyPhen from VEP REST API (default)
    else:    
        # extract the variants from the vcf file
        vcf = get_vcf(args)

        # need to break into chunks of 200 variants - maximum POST size is 200
        def divide_variants_list(list, n):
            for i in range(0, len(list), n):
                yield list[i:i + n]

        variant_list = list(divide_variants_list(vcf, 200)) # returns a list of lists

        # combine all variants in a string to submit to VEP
        print("\nSubmitting variants to VEP server (batch = 200 variants)...")
        n = 0
        for v in variant_list:
            n += 1
            print("processing batch", n)
            variants = ''
            for c in v:
                loc = " ".join(c) 
                variants += "\"" + loc + "\" ," 
            variants = variants[:-1]
            vep_result = vep_output(args, variants, colnames)
            results = pd.concat([results, vep_result])
            time.sleep(2) # avoid taxing server
        print("Done")
    
        # find if mutations are in PDB structures or not
        print("\nGetting PDB IDs...", end='')
        get_pdb_id(results)
        print("Done") 

    # end if args.t

    print("\nGenerating iCn3D URLs...")
    #url_list = defaultdict(dict)
    #gene_id_list = list(set(results.EnsID)) # unique set of gene ids

    #for gene in gene_id_list:
    #    print("=========================\nGetting link for ", gene)
    #    gene_res = results[results['EnsID'] == gene] # just submit the results for a given EnsID
    #    url_list[gene] = get_iCn3D_path(gene_res)
    get_iCn3D_path(results)

    # generate an html page with results dataframe, iCn3D links
    print("\nPrinting html file...")
    print_html(args, results) # url_list, 

    # generate a .csv file with results dataframe, iCn3D links
    print("Printing .csv file...")
    print_csv(args, results) # url_list, 

if __name__ == '__main__':   
    main(cli().parse_args())   
