import os
import sys
import pandas as pd
import obonet
import Bio
from Bio.UniProt.GOA import gafiterator
from datetime import datetime
import ia as ia


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
            # db = entry['DB']
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
    ALL_LISTS = {
        "EXPEC": EXPEC,
        "COMPEC": COMPEC,
        "AUTHEC": AUTHEC,
        "CUREC": CUREC,
        "IEA": IEA
    }
    ALL_CODES = set(code for codes in ALL_LISTS.values() for code in codes)
    input_codes = [code.strip().upper() for code in evidence_codes.split(',')]
    accepted_codes = set()
    for code in input_codes:
        if code in ALL_LISTS:
            accepted_codes.update(ALL_LISTS[code])
        elif code in ALL_CODES:
            accepted_codes.add(code)
        else:
            print(f"Warning: '{code}' is not a recognized evidence code.")
    print(accepted_codes)
    return annotation_df[annotation_df['evidence'].isin(list(accepted_codes))]


def change_date_to_year(annotation_df):
    new_df = annotation_df.copy()
    new_df['year'] = new_df['date'].apply(lambda x: datetime.strptime(x, "%Y%m%d").year)
    del new_df['date']
    return new_df


def divide_by_throughput(df, threshold=100):
    reference_counts = df.groupby('reference')['EntryID'].nunique()
    HTP_reference = reference_counts[reference_counts >= threshold].index
    LTP_reference = reference_counts[reference_counts < threshold].index
    print(len(HTP_reference), len(LTP_reference))
    HTP_data = df[df['reference'].isin(HTP_reference)]
    LTP_data = df[df['reference'].isin(LTP_reference)]        
    return HTP_data.iloc[:,[0,1,2]], LTP_data.iloc[:,[0,1,2]]


def replace_alternate_GO_terms(df, graph):
    # Read the ontology graph (obo file)
    ontology_graph = obonet.read_obo(graph)
    # Create a dictionary to map alternate GO_term back to main GO_term (for example, alt_id: GO:0044822 maps to id: GO:0003723)
    alt_id_map = {id_: data.get('alt_id') for id_, data in ontology_graph.nodes(data=True)}
    alt_id_to_id = {}
    for id_, alt_ids in alt_id_map.items():
        if alt_ids is not None:
            for alt_id in alt_ids:
                alt_id_to_id[alt_id] = id_
    # Replace alternate GO_term with main GO_term
    df['GO_term'] = df['GO_term'].apply(lambda x: alt_id_to_id[x] if x in alt_id_to_id else x)
    return df


def propagate_and_ia(df, graph, year, single_year, throughput_classify=1):
    ontology_graph = ia.clean_ontology_edges(obonet.read_obo(graph))
    roots = {'P': 'GO:0008150', 'C': 'GO:0005575', 'F': 'GO:0003674'}
    subontologies = {aspect: ia.fetch_aspect(ontology_graph, roots[aspect]) for aspect in roots}
    subont_nodes = {aspect: set(subontologies[aspect].nodes()) for aspect in roots}
    
    raw_df = df.loc[:,['UniProt_ID', 'GO_term', 'aspect', 'reference']]
    raw_df.columns = ['EntryID','term','aspect','reference']
    # Note: some leaf terms in the past do not exist in the 2022 graph, remove
    keep = raw_df.apply(lambda row: row['term'] in subont_nodes[row['aspect']], axis=1)
    new_df = raw_df[keep].copy()

    # Propagate terms, new_df has 'reference' column
    print("Propagating terms")
    if throughput_classify > 0:
        annotation_df = ia.propagate_terms(new_df, subontologies)
    else:
        annotation_df = ia.propagate_terms(new_df.iloc[:,[0,1,2]], subontologies)

    # Divide by aspect after propagating
    df = df.loc[:,['UniProt_ID', 'GO_term', 'aspect', 'reference']]
    df.columns = ['EntryID', 'term', 'aspect', 'reference']
    if throughput_classify > 0:
        HTP_annotation_df, LTP_annotation_df = divide_by_throughput(annotation_df)
        HTP_df, LTP_df = divide_by_throughput(new_df)

    # Matching IA with GO terms
    annotation_df = annotation_df.loc[:,['EntryID', 'term', 'aspect']]
    if os.path.exists(f'0_data/ia_{year}.csv'):
        annotation_df.columns = ['UniProt_ID', 'GO_term', 'aspect']
        if throughput_classify > 0:
            HTP_annotation_df.columns = ['UniProt_ID', 'GO_term', 'aspect']
            LTP_annotation_df.columns = ['UniProt_ID', 'GO_term', 'aspect']
    else:
        print('Counting Terms')
        aspect_counts = dict()
        aspect_terms = dict()
        term_idx = dict()
        for aspect, subont in subontologies.items():
            aspect_terms[aspect] = sorted(subont.nodes)  # ensure same order
            term_idx[aspect] = {t:i for i,t in enumerate(aspect_terms[aspect])}
            aspect_counts[aspect] = ia.term_counts(annotation_df[annotation_df.aspect==aspect], term_idx[aspect])
            assert aspect_counts[aspect].sum() == len(annotation_df[annotation_df.aspect==aspect]) + len(aspect_terms[aspect])
        # since we are indexing by column to compute IA, 
        # let's convert to Compressed Sparse Column format
        sp_matrix = {aspect:dok.tocsc() for aspect, dok in aspect_counts.items()}

        # Compute IA
        print('Computing Information Accretion')
        aspect_ia = {aspect: {t:0 for t in aspect_terms[aspect]} for aspect in aspect_terms.keys()}
        for aspect, subontology in subontologies.items():
            for term in aspect_ia[aspect].keys():
                aspect_ia[aspect][term] = ia.calc_ia(term, sp_matrix[aspect], subontology, term_idx[aspect])
        ia_df = pd.concat([pd.DataFrame.from_dict(
            {'term':aspect_ia[aspect].keys(), 
            'ia': aspect_ia[aspect].values(), 
            'aspect': aspect}) for aspect in subontologies.keys()])
        
        # all counts should be non-negative
        assert ia_df['ia'].min() >= 0
        ia_df[['term','ia']].to_csv(f'0_data/ia_{year}.csv', header=True, sep='\t', index=False)

        annotation_df.columns = ['UniProt_ID', 'GO_term', 'aspect']
        
    merge_metrics(new_df.iloc[:,[0,1,2]], annotation_df, year, single_year, "All")
    if throughput_classify > 0:
        merge_metrics(HTP_df, HTP_annotation_df, year, single_year, "HTP")
        merge_metrics(LTP_df, LTP_annotation_df, year, single_year, "LTP")
    return


def get_protein_leaf_metrics(df):
    # For each protein, calculate the number of unique terms it is annotated with in each aspect
    num_unique_terms_aspect = df.groupby(['UniProt_ID', 'aspect'])['GO_term'].unique().apply(lambda x: len(set(x))).reset_index()
    num_unique_terms_aspect.columns = ['UniProt_ID', 'aspect', 'UniCount']
    # For each protein, calculate the number of terms it is annotated with in each aspect
    num_terms_aspect = df.groupby(['UniProt_ID', 'aspect'])['GO_term'].count().reset_index()
    num_terms_aspect.columns = ['UniProt_ID', 'aspect', 'Count']
    num_df = pd.merge(num_terms_aspect, num_unique_terms_aspect, on=['UniProt_ID', 'aspect'])
    return num_df


def get_protein_propagated_metrics(df, year):
    terms_aspect = df.groupby(['UniProt_ID', 'aspect'])['GO_term'].apply(list).reset_index()
    terms_aspect['IA'] = [0 for _ in range(len(terms_aspect))]
    if not os.path.exists(f'0_data/ia_{year}.csv'):
        print("Calculate Information Accretion first")
        return
    else:
        print(f"Using Information Accretion of {year}")
        ia_df = pd.read_csv(f'0_data/ia_{year}.csv', sep="\t", header=0)
        for _,row in terms_aspect.iterrows():
            terms = row['GO_term']
            ia_terms = ia_df[ia_df['term'].isin(terms)]
            ic = ia_terms['ia'].sum()
            terms_aspect.at[_, 'IA'] = ic
    return terms_aspect[['UniProt_ID', 'aspect', 'IA']]
   

def merge_metrics(leaf_df, propagated_df, year, single_year, identifier):
    leaf_df.columns = ['UniProt_ID', 'GO_term', 'aspect']
    leaf_metrics = get_protein_leaf_metrics(leaf_df)
    propagated_metrics = get_protein_propagated_metrics(propagated_df, year)
    merged_metrics = pd.merge(leaf_metrics, propagated_metrics, on=['UniProt_ID', 'aspect'])
    merged_metrics.columns = ['Genes', 'Aspect', 'Count', 'UniCount', 'IA']
    if not os.path.exists("0_data/processed_GO"):
        os.makedirs("0_data/processed_GO")
    merged_metrics.to_csv(f'0_data/processed_GO/{identifier}Proteins_{single_year}.tsv', index=False, sep="\t", header=True)
    return merged_metrics


if __name__ == "__main__":
    if len(sys.argv) != 7:
        print(len(sys.argv))
        print("Usage: script.py <GAF file path> <Year> <Graph file path>")
        sys.exit(1)
    
    single_year = sys.argv[1]
    single_GAF = sys.argv[2]
    pivot_year = int(sys.argv[3])
    pivot_graph = sys.argv[4]
    evidence_codes = sys.argv[5]    
    throughput_classify = sys.argv[6]

    data = extract_GAF_file(single_GAF)
    EXPEC_data = change_date_to_year(filter_evidence_codes(data, evidence_codes))
    EXPEC_data = replace_alternate_GO_terms(EXPEC_data, pivot_graph)
    # EXPEC_data.to_csv(f'0_data/processed_GO/EXPEC_{single_year}.tsv', index=False, sep="\t", header=True)
    filtered_data = EXPEC_data
    propagate_and_ia(filtered_data, graph = pivot_graph, year = int(pivot_year), 
                     single_year=single_year, throughput_classify = int(throughput_classify))
    
    