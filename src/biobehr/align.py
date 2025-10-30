import numpy as np
from pyfaidx import Fasta

# You can do this with AlignIO but it's really slow. Computing alignment.frequencies takes over a minute for large MSA.
# This new implementation with read_msa -> alignment_freqs takes <1s even for large MSA. I checked that the result is the identical.
# bin_nt = (b'-', b'G', b'A', b'T', b'C') # for muscle
bin_nt = (b'-', b'g', b'a', b't', b'c') # for mafft

# Okay and now they're both just single-line lambda functions, after all that. Yeesh. 
read_msa = lambda fa: np.array([np.asarray(r, dtype='|S1') for r in Fasta(fa)])
alignment_freqs          = lambda a   : np.array([np.sum( a == c                   , axis=0) for c in bin_nt])
alignment_freqs_weighted = lambda a, w: np.array([np.sum((a == c).astype(float) * w, axis=0) for c in bin_nt])


# Calculate the entropy at each position in the alignment
def alignment_entropy(freqs):
	probs = freqs / np.sum(freqs[:,0]) # convert frequency counts to probabilities
	masked_probs = np.where(probs > 0, probs, 1e-10) # mask zero probabilities to avoid log(0)
	return -np.sum(probs * np.log2(masked_probs), axis=0) # entropy 


# Note: I tried also smoothing the entropy using sliding windows but then it made the cutoffs for the conserved region less precise. 
def merge_small_high_entropy_segments(below_threshold, max_high_entropy_size):
	# Find the start and end points of high-entropy segments
	high_entropy_segments = np.where(np.diff(np.concatenate(([1], below_threshold, [1]))))[0].reshape(-1, 2)

	# Merge small high-entropy segments
	for start, end in high_entropy_segments:
		if (end - start) <= max_high_entropy_size:
			below_threshold[start:end] = True

	return below_threshold

# max high entropy size is a param that we definitely need to play around with. should maybe only be 2, or even 1
def find_conserved_regions(entropy, threshold=1.0, min_region_size=100, max_high_entropy_size=3):
	# Identify positions below the threshold
	below_threshold = entropy < threshold

	# Merge small high-entropy segments
	below_threshold = merge_small_high_entropy_segments(below_threshold, max_high_entropy_size)

	# Find the start and end points of continuous regions below the threshold
	regions = np.where(np.diff(np.concatenate(([0], below_threshold, [0]))))[0].reshape(-1, 2)

	# Filter regions by minimum size
	valid_regions = regions[(regions[:, 1] - regions[:, 0]) >= min_region_size]

	return valid_regions


str_nt = ( '-',  'G',  'A',  'T',  'C')
map_nt = np.vectorize({i:c for i,c in enumerate(str_nt)}.get)

# Here we define consensus region as the region of high homology (low entropy),
# and the consensus sequence as the most-frequently-observed nt within that region. 
# seq_counts can be a np array of counts for each sequence to weight their frequency by, optionally
def consensus_sequence(msa_file, entropy_threshold=1.0, min_region_size=100, seq_counts=None):
	alignment = read_msa(msa_file)
	if seq_counts is None:
		freqs = alignment_freqs(alignment)  # get the frequency of each base at each position
	else:
		assert isinstance(seq_counts, np.ndarray) and seq_counts.ndim == 1 and seq_counts.shape[0] == alignment.shape[0]
		weights = seq_counts.astype(float)[:,None]
		freqs = alignment_freqs_weighted(alignment, weights)

	# Remove positions whose consensus is gap. We need to do this before we assess regions of low entropy, because
	# rare insertions yield long regions that falely appear to be low-entropy (bc most sequences have gap).
	gaps = np.where(np.argmax(freqs, axis=0) == 0)[0]
	freqs = np.delete(freqs, gaps, axis=1)

	# Figure out the start & end positions of the consensus region
	entropy = alignment_entropy(freqs) # helpful: px.bar(entropy)
	regions = find_conserved_regions(entropy, threshold=entropy_threshold, min_region_size=min_region_size)
	
	# Set start & end to 0,0 (0 length consensus) if no conserved regions are found
	start,end = regions[np.argmax(regions[:,1] - regions[:,0])] if len(regions) > 0 else (0,0)

	# Get the most-frequently-observed nt within the consensus region.
	consensus_idx = np.argmax(freqs, axis=0)[start:end]
	consensus_str = ''.join(map_nt(consensus_idx)) if len(consensus_idx) > 0 else ''

	return dict(
		alignment_length=alignment.shape[1],
		consensus_length=len(consensus_str),
		consensus_start=start,
		consensus_end=end,
		consensus=consensus_str
	)




