import argparse
import requests
import xml.etree.ElementTree as ET
from datetime import datetime
import xmltodict # conda-forge or pip installed
import concurrent.futures
import json
import time
import os
from tqdm import tqdm
from pathlib import Path

# Not the best form to have global vars like this, but it's ok for here & less annoying than
# passing them around through the call stack so many times. 
AUTH_EMAIL = None
API_KEY = None

base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
url = lambda tool: f'{base_url}/{tool}.fcgi'

def params(tool, **kwargs):
	# Required to take in a tool, but anything else is optional. 
	assert AUTH_EMAIL and API_KEY, 'NCBI credentials not initialized'
	p = dict(email=AUTH_EMAIL, api_key=API_KEY, tool=tool, db='sra')
	p.update(kwargs)
	return p

def build_sra_search_query(**kwargs):
	# I don't think the quotes & capitalizations are necessary but w/e
	terms = [f'"{v}"[{k}]' for k,v in kwargs.items()]
	return ' AND '.join(terms)


# =============================================================

def main():
	parser = argparse.ArgumentParser()

	# For now, provide email & key in config yaml file
	parser.add_argument('--email', required=True, help='NCBI user email')
	parser.add_argument('--key', required=True, help='NCBI API Key')	
	parser.add_argument('-o', '--output-dir', help='dir in which to save the three output files.')
	parser.add_argument('-p', '--prefix', default='', help='Name files with a prefix (optional)')
	parser.add_argument('-f', '--force', action='store_true', help='Allow overwrite of existing files.')
	opts, search_terms = parser.parse_known_args()

	# Set NCBI credentials
	global AUTH_EMAIL, API_KEY
	AUTH_EMAIL = opts.email
	API_KEY = opts.key

	# Check/create output directory
	outdir = opts.output_dir if opts.output_dir else f"eutils_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}"
	outdir = Path(outdir).absolute()
	if outdir.exists():
		assert outdir.is_dir(), f'Output dir {outdir} already exists but is not a directory!'
	else:
		os.makedirs(outdir)


	# Check/create output files
	prefix = f'{opts.prefix}_' if opts.prefix and len(opts.prefix) > 0 else ''
	rec_fp = outdir / f'{prefix}records.txt'
	srs_fp = outdir / f'{prefix}samples.jsonl'
	srx_fp = outdir / f'{prefix}experiments.jsonl'
	if not opts.force:
		assert not (rec_fp.exists() or srs_fp.exists() or srx_fp.exists()), 'Output file(s) already exist!'

	# Parse remaining inputs as search terms & build query
	search_terms = parse_search_terms(search_terms)
	query = build_sra_search_query(**search_terms)
	divider = '\n' + '-'*18 + '\n'
	print(f'{divider}Entrez Query Plan: {query}\nOutput dir: {outdir}{divider}')

	# Execute queries
	num_results = get_num_results(query)
	record_ids = esearch_get_record_ids(query, num_results, rec_fp)
	parallel_efetch_records(record_ids, srs_fp, srx_fp)

	# Why separate into two ndjson files rather than one?
	# It's not strictly necessary, but there will be fewer unique samples than experiment/runs.


def parse_search_terms(terms):
	assert len(terms) % 2 == 0, f'Invalid search terms {terms}'
	return {
		terms[i].lstrip('--').title():terms[i+1] 
		for i in range(0, len(terms), 2)
	}

# E.g. if you already have a records.txt file and want to efetch on all or some batches
def parse_record_ids(record_ids_fp):
	with open(record_ids_fp) as f:
		return [x for x in (l.strip() for l in f) if len(x) > 0 and not x.startswith('#')]

'''
efm_long_read = build_sra_search(
	Organism = 'Enterococcus faecium',
	Platform = 'oxford nanopore',
	Source = 'genomic',
	Strategy = 'wgs'
)

staph_short_read = build_sra_search(
	Organism = 'Staphylococcus aureus',
	Platform = 'illumina',
	Source = 'genomic',
	Strategy = 'wgs'
)

staph_epi_short_read = build_sra_search(
	Organism = 'Staphylococcus epidermidis',
	Platform = 'illumina',
	Source = 'genomic',
	Strategy = 'wgs'
)
'''

### Figure out how many results there will be
def get_num_results(query):
	print('Get esearch query scope')
	search_params = params('esearch', rettype='count', retmode='json', term=query)
	r = requests.get(url('esearch'), params=search_params)
	return int(r.json()['esearchresult']['count']) # 152832

### Obtain all the record IDs with esearch (up to 10k at a time, which is the limit):
def esearch_get_record_ids(query, num_results, output_file=None):
	print(f'Fetch {num_results} record IDs with esearch')
	record_ids = []
	step = 10000 # 10k is the esearch query limit, so we can do 10k at a time

	with tqdm(total=num_results, desc='ESearch', unit='records') as pbar:
		for i in range(0, num_results, step):
			search_params = params(
				tool = 'esearch',
				retmode = 'json',
				term = query,
				retstart = i,
				retmax = step
			)
			r = requests.get(url('esearch'), params=search_params)
			assert r.status_code == 200, r.text # print response if error status
			batch_ids = r.json()['esearchresult']['idlist']
			record_ids += batch_ids
			pbar.update(len(batch_ids))
	
	assert len(record_ids) == len(set(record_ids)), 'Duplicate records'
	assert len(record_ids) == num_results, 'Missing records'

	if output_file:
		print(f'Save record IDs to file at {output_file}')
		with open(output_file, 'w') as f:
			f.write(f'# ESearch run at {datetime.now()}\n')
			f.write(f'# Query: {query}\n')
			f.write(f'# Results: {num_results}\n')
			f.write('\n'.join(record_ids))

	return record_ids


### Obtain all the data info with efetch (1k at a time or in parallel, it's pretty slow)

# efetch can't do json retmode, so you need to parse the XML, unfortunately. 

# fetch_data = {'id': r.json()['esearchresult']['idlist']}
# r = requests.post(url('efetch'), params=params('efetch'), data={'id': record_ids})
# root = ET.fromstring(r.content)

# NOTE: this will need to be modified depending on what info you care about & possibly 
# the format of the data. I've tested it with several different types of query results but YMMV.
def parse_efetch_result(experiment_package):
	s = experiment_package.find('SAMPLE')
	x = experiment_package.find('EXPERIMENT')
	r = experiment_package.find('RUN_SET')

	# Sample
	biosample = s.find('./IDENTIFIERS/EXTERNAL_ID[@namespace="BioSample"]')
	scientific_name = s.find('./SAMPLE_NAME/SCIENTIFIC_NAME')
	
	sample = dict(
		SRS = s.attrib.get('accession'),
		biosample = biosample.text if isinstance(biosample, ET.Element) else '',
		scientific_name = scientific_name.text if isinstance(scientific_name, ET.Element) else ''
	)
	
	# Very rarely, no attributes. Usually these are invalid anyway. 
	if isinstance(s.find('./SAMPLE_ATTRIBUTES'), ET.Element):
		attrs = xmltodict.parse(ET.tostring(s.find('./SAMPLE_ATTRIBUTES')))['SAMPLE_ATTRIBUTES'].get('SAMPLE_ATTRIBUTE')
		if isinstance(attrs, dict): attrs = [attrs] # edge case of 1 attr
		for attr in attrs: sample[attr['TAG']] = attr.get('VALUE') # weirdly there are a few that have a tag but no value, so we'll just make those null. 

	# Extract platform info in advance, so it can be used to determine instrument model below.
	platform = [child.tag for child in x.find('PLATFORM')]
	assert len(platform) == 1, f'Found multiple platforms in single experiment package: {platform}'

	# Experiment & Runs
	experiment_run = dict(
		SRX = x.attrib.get('accession'),
		title = x.find('./TITLE').text,
		SRS = x.find('./DESIGN/SAMPLE_DESCRIPTOR').attrib.get('accession'),
		# This would be nice to do to see the bioproject diversity across samples, but need to look into it more bc it's not consistent across all data. 
		# bioproject = x.find('./STUDY_REF/IDENTIFIERS/EXTERNAL_ID[@namespace="BioProject"]').text,
		library_layout = x.find('./DESIGN/LIBRARY_DESCRIPTOR/LIBRARY_LAYOUT/*').tag,
		design_description = x.find('./DESIGN/DESIGN_DESCRIPTION').text,

		# instrument_model = x.find('PLATFORM/ILLUMINA/INSTRUMENT_MODEL').text,
		# I looked at the oputput of this and model is the only thing, so this doesn't provide more info. 
		# But that's only true depending on the platform. Modified to be a bit more forgiving:
		platform = platform[0],
		instrument_model = x.find(f'PLATFORM/{platform[0]}/INSTRUMENT_MODEL').text,

		# run information
		runs = int(r.attrib.get('runs')),
		bases = int(r.attrib.get('bases')),
	
		# Note there can be multiple runs for a given Experiment; this can be just like due to 
		# filesize or whatever. You just combine them before assembly. I think also you can from sra
		# via the experiment accession rather than the run accession. 
		SRR = ','.join([x.attrib.get('accession') for x in r.findall('./RUN')])
	)

	return sample, experiment_run


def efetch(ids):
	r = requests.post(url('efetch'), params=params('efetch'), data={'id': ids})
	assert r.status_code == 200, r.content
	return [parse_efetch_result(e) for e in ET.fromstring(r.content)]


### Efetch (Process in parallel)
def parallel_efetch_records(record_ids, srs_fp, srx_fp):
	chunks = [record_ids[i:i+1000] for i in range(0, len(record_ids), 1000)]

	# API allows 10 requests per second, we'll limit to 8 at a time but they'll take longer than a second anyway
	with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
		
		print('Submit rate-limited jobs')
		futures = {}
		for i,chunk in enumerate(tqdm(chunks, total=len(chunks), desc='Submit', unit='chunks')):
			futures[executor.submit(efetch, chunk)] = i
			if i % 10 == 0: time.sleep(1.2) # delay >1 second every 10 requests so we don't overwhelm the API

		print('Waiting for jobs to complete (this can take a while to update the first time)...')
		with tqdm(total=len(chunks), desc='Complete', unit='chunks') as pbar:
			for future in concurrent.futures.as_completed(futures):
				chunk_idx = futures[future]
				try:
					with open(srs_fp, 'a') as s:
						with open(srx_fp, 'a') as e:
							for sample, experiment_run in future.result():
								s.write(json.dumps(sample) + '\n')
								e.write(json.dumps(experiment_run) + '\n')
					pbar.update(1)
				except Exception as e:
					# Report error w/o breaking the pbar
					tqdm.write(f'Chunk {chunk_idx} generated an exception: {e}')

'''
When I ran this for 9.5k records on M4:
- job submission was almost instant
- took 1-2 minutes (waiting for jobs to complete) before the 1st job completed
- took another minute as progress went through to complete all the jobs
'''

# Jobs can sometimes fail for various reasons, but you may not want to start 
# from scratch. This way you can manually rerun a few chunks that failed 
# (this code is one-off based on the Sa run that I did, but keeping just for 
# records, in case you want to do something like this again). For example:
# error_chunks = [52, 139, 140, 150, 151]

def rerun_error_chunks(chunks, error_chunk_idxs):
	for chunk_id in error_chunk_idxs:
		print(f'Attempt chunk {chunk_id}')
		result = efetch(chunks[chunk_id])
		with open('sample_chunks_rerun.jsonl', 'a') as s:
			with open('experiment_chunks_rerun.jsonl', 'a') as e:
				for sample, experiment_run in result:
					s.write(json.dumps(sample) + '\n')
					e.write(json.dumps(experiment_run) + '\n')
		print(f'Succeed chunk {chunk_id}')



if __name__ == '__main__':
	main()