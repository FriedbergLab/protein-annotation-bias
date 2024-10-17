import os
import sys
import pandas as pd
import Bio
from Bio.UniProt.GOA import gafiterator
from datetime import datetime
import subprocess

EXPEC = [
    "EXP",
    "IDA",
    "IPI",
    "IMP",
    "IGI",
    "IEP",
    "HTP",
    "HDA",
    "HMP",
    "HGI",
    "HEP"
]

COMPEC = [
    "ISS",
    "ISO",
    "ISA",
    "ISM",
    "IGC",
    "IBA",
    "IBD",
    "IKR",
    "IRD",
    "RCA"
]

AUTHEC = [
    "TAS",
    "NAS"
]

CUREC = [
    "IC",
    "ND"
]

IEA = ["IEA"]


def process_gaf_file(gaf_file):
    with open(gaf_file, 'r') as f:
        content = f.read()
    
    # Find the position of "!gaf-version" line
    gaf_version_pos = content.find("!gaf-version")
    
    if gaf_version_pos == -1:
        # If "!gaf-version" is not found, process from the beginning
        return content
    else:
        # If "!gaf-version" is found, process from that line onwards
        return content[gaf_version_pos:]
    

def extract_GAF_file(gaf_file_path):
    data = []
    # Open the GAF file and parse it
    gaf_content = process_gaf_file(gaf_file_path)
    temp_gaf_file = gaf_file_path + '.temp'
    with open(temp_gaf_file, 'w') as f:
        f.write(gaf_content)
    with open(temp_gaf_file, 'r') as temp_file:
        # Use gafiterator to parse the GAF file
        for entry in gafiterator(temp_file):
            if 'NOT' in entry['Qualifier']:
                continue
            # Extract the required columns
            UniProt_ID = entry['DB_Object_ID']
            UniProt_GeneName = entry['DB_Object_Symbol']
            GO_term = entry['GO_ID']
            aspect = entry['Aspect']
            reference = entry['DB:Reference'][0]
            evidence = entry['Evidence']
            date = entry['Date']
            # Store the extracted information
            data.append([UniProt_ID, UniProt_GeneName, GO_term, aspect, reference, evidence, date])
    os.remove(temp_gaf_file)
    return pd.DataFrame(data, columns=['UniProt_ID', 'UniProt_GeneName', 'GO_term', 'aspect', 'reference','evidence', 'date'])


def filter_evidence_codes(annotation_df, evidence_codes):
    return annotation_df[annotation_df['evidence'].isin(evidence_codes)]


def change_date_to_year(annotation_df):
    new_df = annotation_df.copy()
    new_df['year'] = new_df['date'].apply(lambda x: datetime.strptime(x, "%Y%m%d").year)
    del new_df['date']
    return new_df

def EDirect_get_year(ref):
    try:
        # Run the efetch and xtract commands
        result = subprocess.run(
            ['efetch', '-db', 'pubmed', '-id', str(ref), '-format', 'docsum'],
            stdout=subprocess.PIPE, text=True
        )
        # Extract the PubDate using xtract
        output = subprocess.run(
            ['xtract', '-pattern', 'DocumentSummary', '-element', 'PubDate'],
            input=result.stdout, stdout=subprocess.PIPE, text=True
        )
        # Process the output to extract the year
        pub_date = output.stdout.strip()
        if pub_date:
            year = pub_date.split()[0]  # Extract the year if multiple fields are returned
            if year.isdigit():  # Check if the extracted string is a valid year
                return year
            else:
                return '0'
        else:
            return '0'
    except Exception as e:
        print(f"Error fetching publication year for PMID {ref}: {e}")
        return '0'

def get_pubyear(ref_data, checkpoint_threshold, endfile):
    checkpoint = endfile.replace(".csv","_checkpoint.csv")
    checkpoint_interval = int(checkpoint_threshold)
    checkpoint_index = 0
    if os.path.exists(checkpoint):
        PMID_year_df = pd.read_csv(checkpoint, sep="\t",header=0)
        PMID_year = dict(zip(PMID_year_df["PMID"].astype(str), PMID_year_df["pub_year"].astype(str)))
        processed = set(PMID_year.keys())    
    else:
        PMID_year = dict()
        processed = set()
    
    for ref in ref_data:
        if 'PMID' in ref:
            ID = str(ref.split(":")[-1])
            if ID not in processed:            
                pub_year = EDirect_get_year(ID)
                PMID_year[ID]=pub_year
                processed.add(ID)
                checkpoint_index+=1
            if checkpoint_index > 0 and checkpoint_index % checkpoint_interval == 0:
                checkpoint_df = pd.DataFrame(PMID_year.items(), columns=["PMID","pub_year"])
                checkpoint_df.to_csv(checkpoint,sep="\t",index=False,header=True)
                print(f"Checkpoint saved after {checkpoint_index} entries")
    
    PMID_year_df = pd.DataFrame(PMID_year.items(),columns=["PMID","pub_year"])
    PMID_year_df.to_csv(endfile,sep="\t",header=True,index=False)
    return PMID_year_df


def merge_pubyear(full_data, PMID_year_df, df_provided=True):
    if df_provided:
        PMID_year = PMID_year_df
    else:   
        if os.path.exists(PMID_year_df):
            PMID_year = pd.read_csv(PMID_year_df,sep="\t",header=0,dtype={"PMID":str,"pub_year":str})
        else:
            print("Find publication years for PMID first")
    full_data['PMID']=full_data['reference'].apply(lambda x: str(x.split(":")[-1]) if "PMID" in x else x)
    merge_df = pd.merge(full_data,PMID_year,on=["PMID"],how="left")
    print(merge_df.head())
    merge_df.to_csv(f'3_delay_curation/annotation_pub_years.csv',sep="\t",header=True,index=False)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 3.1_extract_pubyear.py <GAF_file>")
        sys.exit(1)
    single_GAF = sys.argv[1]
    
    data = extract_GAF_file(single_GAF)
    EXPEC_data = change_date_to_year(filter_evidence_codes(data, EXPEC))
    ref_data = EXPEC_data["reference"].unique().tolist()
    print(f"Total: {len(ref_data)} unique PMID articles")
    PMID_year_df = get_pubyear(ref_data, 100, endfile="3_delay_curation/PMID_year_EDirect.csv")
    merge_pubyear(EXPEC_data, PMID_year_df, df_provided=True)
    # Use this if the PMID_year_df is already saved in a file:
    # merge_pubyear(EXPEC_data, "3_delay_curation/PMID_year_EDirect.csv", df_provided=False) 
    


