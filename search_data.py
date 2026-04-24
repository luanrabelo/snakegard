import argparse
import asyncio
import logging
import time
import os

import requests

import pandas as pd

from Bio import Entrez
from dotenv import load_dotenv

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y/%m/%d - %H:%M:%S',
)

def load_api_keys(**kwargs):
    env_file = kwargs.get('env_file', 'apikeys.env')
    verbose = kwargs.get('verbose', True)

    if os.path.exists(env_file):
        if verbose:
            logging.info(f"Loading API keys from {env_file}")
        load_dotenv(env_file)
    else:
        if verbose:
            logging.warning(f"API keys file '{env_file}' not found, checking environment variables")
        load_dotenv()

    api_keys = {
        'IUCN_API_KEY': os.getenv('IUCN_API_KEY'),
    }
    
    return api_keys

def search_organelle_genomes(**kwargs):
    species = kwargs.get("species_name", None)
    email = kwargs.get("email", None)
    verbose = kwargs.get("verbose", False)
    organelle = kwargs.get("organelle", "mitochondrion")

    if not species:
        if verbose:
            logging.error('Species name is required for searching the NCBI Nucleotide database.')
        return None

    if email:
        Entrez.email = email

    label = "Mitochondrion" if organelle == "mitochondrion" else "Chloroplast"

    if verbose:
        logging.info(f'Searching for {species} in NCBI Nucleotide ({label})...')

    query = f'("{species}"[Organism]) AND ("genome"[Title]) AND ("{organelle}"[Title]) AND ("complete"[Title] OR "partial"[Title]) AND ({organelle}[filter])'

    handle = Entrez.esearch(db='nucleotide', term=query, retmax=100000)
    record = Entrez.read(handle)
    result = str(record['Count'])
    handle.close()

    if verbose:
        logging.info(f'Found {result} results for {species} in NCBI Nucleotide ({label}).')

    return result

def search_sra(**kwargs):
    species = kwargs.get("species_name", None)
    email = kwargs.get("email", None)
    verbose = kwargs.get("verbose", False)
    read_type = kwargs.get("read_type", "short")

    if not species:
        if verbose:
            logging.error('Species name is required for searching the NCBI SRA database.')
        return None

    if email:
        Entrez.email = email

    if read_type == "long":
        label = "SRA Long Reads"
        query = f'("{species}"[Organism]) AND ("biomol rna"[Properties] AND ("platform oxford nanopore"[Properties] OR "platform pacbio smrt"[Properties]))'
    else:
        label = "SRA Short Reads"
        query = f'("{species}"[Organism]) AND ("biomol rna"[Properties] AND "library layout paired"[Properties] AND "platform illumina"[Properties])'

    if verbose:
        logging.info(f'Searching for {species} in NCBI {label}...')

    handle = Entrez.esearch(db='sra', term=query, retmax=100000)
    record = Entrez.read(handle)
    result = str(record['Count'])
    handle.close()

    if verbose:
        logging.info(f'Found {result} results for {species} in NCBI {label}.')

    return result

def get_iucn_data(**kwargs):
    species = kwargs.get('species_name', "Rhinella marina")
    token = kwargs.get('token', None)
    verbose = kwargs.get('verbose', False)
    max_retries = kwargs.get('max_retries', 5)
    retry_delay = kwargs.get('retry_delay', 5)

    if not species:
        if verbose:
            logging.error('Species name is required for searching the IUCN Red List API.')
        return None
    
    if not token:
        if verbose:
            logging.warning('IUCN API token not provided. You may encounter rate limits or failed requests.')
        return None
    
    if verbose:
        logging.info(f'Searching for {species} in IUCN Red List API...')
    
    data_dict = {
        'Kingdom': '-', 'Phylum': '-', 'Class': '-', 'Order': '-', 
        'Family': '-', 'Genus': str(species).split()[0], 'Species Name': species,
        'Authority': '-', 'Main Common Name': '-', 'Red List Category': '-',
        'Population Trend': '-', 'Threats List': '-', 'Possibly Extinct': 'No',
        'Possibly Extinct In Wild': 'No'
    }

    species_parts = species.split()
    if len(species_parts) < 2:
        if verbose: 
            logging.error(f"IUCN - {species}: Species name must include at least genus and species epithet.")
        return data_dict

    genus_name, species_epithet = species_parts[0], species_parts[1]
    url = f'https://api.iucnredlist.org/api/v4/taxa/scientific_name?genus_name={genus_name}&species_name={species_epithet}'
    headers = {'accept': 'application/json'}
    if token:
        headers['Authorization'] = token

    for attempt in range(max_retries):
        try:
            response = requests.get(url, headers=headers, timeout=10)
            
            if response.status_code == 200:
                data = response.json()
                taxon_info = data.get('taxon')

                if taxon_info:
                    data_dict['Kingdom'] = str(taxon_info.get('kingdom_name', '-')).capitalize()
                    data_dict['Phylum'] = str(taxon_info.get('phylum_name', '-')).capitalize()
                    data_dict['Class'] = str(taxon_info.get('class_name', '-')).capitalize()
                    data_dict['Order'] = str(taxon_info.get('order_name', '-')).capitalize()
                    data_dict['Family'] = str(taxon_info.get('family_name', '-')).capitalize()
                    data_dict['Genus'] = str(taxon_info.get('genus_name', genus_name)).capitalize()
                    data_dict['Authority'] = taxon_info.get('authority', '-')

                    cn = taxon_info.get('common_names', [])
                    for name_data in cn:
                        if str(name_data.get('main', '')).lower() == 'true':
                            lang = name_data.get('language', '')
                            data_dict['Main Common Name'] = f"{name_data.get('name')} ({lang})" if lang else name_data.get('name')
                            break

                    assessments = data.get('assessments', [])
                    if assessments:
                        global_as = [a for a in assessments if any('global' in str(s.get('description', {}).get('en', '')).lower() for s in a.get('scopes', []))]
                        latest = max(global_as or assessments, key=lambda x: int(x.get('year_published', 0)))
                        
                        codes = {'EX': 'Extinct', 'EW': 'Extinct in the Wild', 'CR': 'Critically Endangered', 'EN': 'Endangered', 'VU': 'Vulnerable', 'NT': 'Near Threatened', 'LC': 'Least Concern', 'DD': 'Data Deficient'}
                        c_code = latest.get('red_list_category_code', '-')
                        data_dict['Red List Category'] = codes.get(c_code, c_code)
                        data_dict['Possibly Extinct'] = 'Yes' if latest.get('possibly_extinct') else 'No'

                        a_id = latest.get('assessment_id')
                        if a_id:
                            detail_url = f'https://api.iucnredlist.org/api/v4/assessment/{a_id}'
                            detail_resp = requests.get(detail_url, headers=headers, timeout=10)
                            if detail_resp.status_code == 200:
                                a_data = detail_resp.json()
                                
                                trend = a_data.get('population_trend')
                                if trend:
                                    t_desc = trend.get('description', {}).get('en', '-') if isinstance(trend.get('description'), dict) else trend.get('description', '-')
                                    data_dict['Population Trend'] = str(t_desc).capitalize()
                                
                                threats = [f"{t.get('description', {}).get('en', '')}" for t in a_data.get('threats', []) if t.get('code')]
                                data_dict['Threats List'] = '; '.join(threats) if threats else '-'

                    if verbose: 
                        logging.info(f"Successfully retrieved IUCN data for {species}.")
                    return data_dict

            elif response.status_code == 429:
                time.sleep(retry_delay * 2)
            else:
                if verbose: 
                    logging.warning(f"Error {response.status_code} for {species} in IUCN API")

        except Exception as e:
            if verbose: 
                logging.error(f"Connection error for {species} in IUCN API: {e}")
        
        time.sleep(retry_delay)

    if verbose: 
        logging.warning(f"Failed to retrieve IUCN data for {species} after {max_retries} attempts.")
    return data_dict


class NCBIRateLimiter:
    """Serializes NCBI calls and enforces ~350ms spacing (max 3 req/s)."""
    def __init__(self):
        self._lock = asyncio.Lock()

    async def call(self, func, **kwargs):
        async with self._lock:
            result = await asyncio.to_thread(func, **kwargs)
            await asyncio.sleep(0.5)
            return result


async def process_species(species_name, ncbi_limiter, iucn_token, verbose, counter, total_pending, already_processed, total_species):
    """Process a single species: 4 NCBI searches + 1 IUCN search concurrently."""
    current = counter['value']
    counter['value'] += 1

    mito_task = ncbi_limiter.call(search_organelle_genomes, species_name=species_name, organelle="mitochondrion", verbose=False)
    chloro_task = ncbi_limiter.call(search_organelle_genomes, species_name=species_name, organelle="chloroplast", verbose=False)
    sra_short_task = ncbi_limiter.call(search_sra, species_name=species_name, read_type="short", verbose=False)
    sra_long_task = ncbi_limiter.call(search_sra, species_name=species_name, read_type="long", verbose=False)
    iucn_task = asyncio.to_thread(get_iucn_data, species_name=species_name, token=iucn_token, verbose=False)

    genome_mito, genome_chloro, sra_short, sra_long, iucn_data = await asyncio.gather(
        mito_task, chloro_task, sra_short_task, sra_long_task, iucn_task
    )

    mito_val = genome_mito if genome_mito is not None else '0'
    chloro_val = genome_chloro if genome_chloro is not None else '0'
    short_val = sra_short if sra_short is not None else '0'
    long_val = sra_long if sra_long is not None else '0'
    iucn_status = iucn_data.get('Red List Category', '-') if iucn_data else 'Error'

    logging.info(
        f'[{current}/{total_pending}] ({already_processed + current}/{total_species}) Processing: {species_name}\n'
        f'\tNCBI Mitochondrion: {mito_val}\n'
        f'\tNCBI Chloroplast:   {chloro_val}\n'
        f'\tSRA Short Reads:    {short_val}\n'
        f'\tSRA Long Reads:     {long_val}\n'
        f'\tIUCN Red List:      {iucn_status}\n'
    )

    if iucn_data is None:
        iucn_data = {
            'Kingdom': '-', 'Phylum': '-', 'Class': '-', 'Order': '-',
            'Family': '-', 'Genus': str(species_name).split()[0], 'Species Name': species_name,
            'Authority': '-', 'Main Common Name': '-', 'Red List Category': '-',
            'Population Trend': '-', 'Threats List': '-', 'Possibly Extinct': 'No',
            'Possibly Extinct In Wild': 'No'
        }

    return [
        iucn_data.get('Kingdom', '-'),
        iucn_data.get('Phylum', '-'),
        iucn_data.get('Class', '-'),
        iucn_data.get('Order', '-'),
        iucn_data.get('Family', '-'),
        iucn_data.get('Genus', '-'),
        iucn_data.get('Species Name', species_name),
        iucn_data.get('Authority', '-'),
        iucn_data.get('Main Common Name', '-'),
        iucn_data.get('Red List Category', '-'),
        iucn_data.get('Population Trend', '-'),
        iucn_data.get('Threats List', '-'),
        iucn_data.get('Possibly Extinct', '-'),
        iucn_data.get('Possibly Extinct In Wild', '-'),
        mito_val,
        chloro_val,
        short_val,
        long_val,
    ]


async def run(input_file, output_file, email, env_file, verbose):
    # 1. Carregar chaves e configurar Entrez
    api_keys = load_api_keys(env_file=env_file, verbose=verbose)
    if email:
        Entrez.email = email
    else:
        logging.warning("NCBI email not provided. Queries might be limited or blocked.")

    # 2. Ler arquivo de entrada
    if not os.path.exists(input_file):
        logging.error(f"Input file '{input_file}' not found.")
        return

    df_input = pd.read_csv(input_file, header=None, names=['species_names'])
    species_list = df_input['species_names'].str.strip().sort_values().unique()

    # 3. Checkpointing
    processed_species = set()

    if os.path.exists(output_file):
        try:
            existing_df = pd.read_csv(output_file, sep='\t', usecols=['Species Name'])
            processed_species = set(existing_df['Species Name'].astype(str).tolist())
            if verbose:
                logging.info(f"Resuming: {len(processed_species)} species already processed.")
        except Exception as e:
            logging.warning(f"Could not parse existing output file ({e}). Starting fresh.")
    else:
        output_dir = os.path.dirname(output_file)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True, mode=0o755)

        header = [
            'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus',
            'Species Name', 'Authority', 'Main Common Name', 'Red List Category',
            'Population Trend', 'Threats List', 'Possibly Extinct',
            'Possibly Extinct In Wild', 'Mitochondrion Count', 'Chloroplast Count',
            'SRA Short Reads Count', 'SRA Long Reads Count'
        ]
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write('\t'.join(header) + '\n')

    pending_species = [s for s in species_list if s not in processed_species]
    total_pending = len(pending_species)
    total_species = len(species_list)

    if verbose:
        logging.info(f"Total species: {total_species} | Already processed: {len(processed_species)} | Pending: {total_pending}")

    if total_pending == 0:
        logging.info("Nothing to process.")
        return

    # 4. Processamento assíncrono
    ncbi_limiter = NCBIRateLimiter()
    species_semaphore = asyncio.Semaphore(10)
    file_lock = asyncio.Lock()
    counter = {'value': 1}

    async def process_and_write(species_name):
        async with species_semaphore:
            try:
                row_data = await process_species(
                    species_name, ncbi_limiter,
                    api_keys.get('IUCN_API_KEY'), verbose,
                    counter, total_pending, len(processed_species), total_species
                )
                sanitized_row = [str(item).replace('\t', ' ') for item in row_data]
                async with file_lock:
                    with open(output_file, 'a', encoding='utf-8') as tsv:
                        tsv.write('\t'.join(sanitized_row) + '\n')
            except Exception as e:
                logging.error(f"Critical error processing {species_name}: {e}")

    tasks = [process_and_write(sp) for sp in pending_species]
    await asyncio.gather(*tasks)

    # 5. Ordenar o arquivo de saída por taxonomia
    logging.info("Sorting output file by taxonomic classification...")
    df_output = pd.read_csv(output_file, sep='\t')
    df_output.sort_values(
        by=['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species Name'],
        inplace=True
    )
    df_output.to_csv(output_file, sep='\t', index=False)

    logging.info("Process completed successfully.")


if __name__ == "__main__":
    argparser = argparse.ArgumentParser(
        description='Search for Organelle Genomes, RNA-seq, and Conservation Data in NCBI databases/IUCN Red List.'
    )

    argparser.add_argument(
        '-i', '--input',
        default='species_brazil.txt',
        type=str,
        help='Input file containing species names, one per line (required).'
    )

    argparser.add_argument(
        '-o', '--output',
        type=str,
        default='search_results.tsv',
        help='Output TSV file to save the search results (default: search_results.tsv).'
    )

    argparser.add_argument(
        '-e', '--email',
        type=str,
        help='Email address to use with NCBI Entrez API (required for large queries).'
    )

    argparser.add_argument(
        '-env',
        '--env_file',
        type=str,
        default="apikeys.env",
        help="Path to the .env file with IUCN_API_KEY."
    )

    argparser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose logging.'
    )

    input_args = argparser.parse_args()

    if os.path.isdir(input_args.output):
        logging.error(f"Output path '{input_args.output}' is a directory. Please provide a file path (e.g., ./search_results.tsv).")
        exit(1)

    asyncio.run(run(
        input_file=input_args.input,
        output_file=input_args.output,
        email=input_args.email,
        env_file=input_args.env_file,
        verbose=input_args.verbose,
    ))