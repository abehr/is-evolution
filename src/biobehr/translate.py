import numpy as np

# Precompute the amino acid array for fast lookup
# Map codons to indices by encoding nucleotides as numbers
_base_to_num = np.zeros(256, dtype=np.uint8)
_base_to_num[ord('A')] = 0
_base_to_num[ord('C')] = 1
_base_to_num[ord('G')] = 2
_base_to_num[ord('T')] = 3

# Initialize the amino acid array with 'X' for unknown codons
aa_array = np.array(['X'] * 64)

# Codon to amino acid mapping
codon_table = {
	'TTT': 'F', 'TTC': 'F',
	'TTA': 'L', 'TTG': 'L',
	'CTT': 'L', 'CTC': 'L',
	'CTA': 'L', 'CTG': 'L',
	'ATT': 'I', 'ATC': 'I',
	'ATA': 'I',
	'ATG': 'M',
	'GTT': 'V', 'GTC': 'V',
	'GTA': 'V', 'GTG': 'V',
	'TCT': 'S', 'TCC': 'S',
	'TCA': 'S', 'TCG': 'S',
	'AGT': 'S', 'AGC': 'S',
	'CCT': 'P', 'CCC': 'P',
	'CCA': 'P', 'CCG': 'P',
	'ACT': 'T', 'ACC': 'T',
	'ACA': 'T', 'ACG': 'T',
	'GCT': 'A', 'GCC': 'A',
	'GCA': 'A', 'GCG': 'A',
	'TAT': 'Y', 'TAC': 'Y',
	'TAA': '*', 'TAG': '*', 'TGA': '*',  # Stop codons
	'CAT': 'H', 'CAC': 'H',
	'CAA': 'Q', 'CAG': 'Q',
	'AAT': 'N', 'AAC': 'N',
	'AAA': 'K', 'AAG': 'K',
	'GAT': 'D', 'GAC': 'D',
	'GAA': 'E', 'GAG': 'E',
	'TGT': 'C', 'TGC': 'C',
	'TGG': 'W',
	'CGT': 'R', 'CGC': 'R',
	'CGA': 'R', 'CGG': 'R',
	'AGA': 'R', 'AGG': 'R',
	'GGT': 'G', 'GGC': 'G',
	'GGA': 'G', 'GGG': 'G',
}

# Precompute the amino acid array for all codons
for codon, aa in codon_table.items():
	idx = (
		(_base_to_num[ord(codon[0])] << 4) +
		(_base_to_num[ord(codon[1])] << 2) +
		_base_to_num[ord(codon[2])]
	)
	aa_array[idx] = aa

# Translate nt -> aa string. np vectorized for high performance.
def translate(seq, to_string=True, continue_after_stop=False):
	
	# Trim sequence to a multiple of 3 and convert to str
	seq = seq[:(len(seq) // 3) * 3]

	# Convert sequence to np array of nt indices
	seq_arr = np.frombuffer(seq.encode('ascii'), dtype=np.uint8)
	seq_nums = _base_to_num[seq_arr]

	# Calculate codon indices & translate to amino acids
	codon_indices = ( (seq_nums[0::3] << 4) + (seq_nums[1::3] << 2) + seq_nums[2::3] )
	aa = aa_array[codon_indices]

	# Optionally, trim to the first stop codon
	if not continue_after_stop:
		stop_idx = np.where(aa == '*')[0]
		if stop_idx.size > 0: aa = aa[:stop_idx[0]+1]

	return ''.join(aa) if to_string else aa


'''
If you store translated orfs as np char array -> tobytes(), this helper function converts back to str.
For example:
orf_bytes = translate(seq, to_string=False).tobytes()
orf_str = from_bytes(orf_bytes)
'''
bytes_to_str = lambda x: ''.join(chr(c) for c in np.frombuffer(x, dtype='<u4'))