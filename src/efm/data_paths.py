from __future__ import annotations
from dataclasses import dataclass, field, fields
from pathlib import Path


# Helper function: declare Path field from basedir / parts
rel = lambda *pathparts: field(init=False, metadata=dict(kind='relpath', data=pathparts))

# Helper function: declare nested dataclass section (w/ same basedir as parent)
section = lambda cls: field(init=False, metadata=dict(kind='section', data=cls))

# Post-init path resolver for both rel & section. This allows them to wait 
# to be resolved until after they get a basedir/root passed in
class AutoResolve:
	def _resolve(self):
		basedir = getattr(self, 'basedir')
		for f in fields(self):
			# Attrs with no metadata will be skipped (but need to handle them correctly)
			meta = f.metadata or {}
			kind = meta.get('kind')
			data = meta.get('data')
			if kind == 'relpath':
				object.__setattr__(self, f.name, basedir.joinpath(*data))
			elif kind == 'section':
				subcls = data
				object.__setattr__(self, f.name, subcls(basedir))

# -----------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class Ncbi(AutoResolve):
	basedir: Path
	complete_genomes: Path = rel('ncbi_complete_datasets.csv')
	nanopore_samples: Path = rel('ncbi_nanopore_samples.csv')
	combined_qc_filtered: Path = rel('ncbi_all_filtered_samples.csv')
	combined_final_samples: Path = rel('ncbi_all_filtered_samples_full.csv')
	def __post_init__(self): self._resolve()

@dataclass(frozen=True, slots=True)
class Shc(AutoResolve):
	basedir: Path
	iso: Path = rel('shc_samples_qc.csv')
	other_ente: Path = rel('other_ente_samples.csv')
	lr_meta: Path = rel('lr_meta_samples.csv')
	sequencing_stats: Path = rel('plasmidsaurus_sequencing_stats.csv')
	def __post_init__(self): self._resolve()

@dataclass(frozen=True, slots=True)
class IseScan(AutoResolve):
	basedir: Path
	summary: Path = rel('isescan_summary.csv')
	per_contig: Path = rel('isescan_per_contig.csv')
	def __post_init__(self): self._resolve()

@dataclass(frozen=True, slots=True)
class LrMetaSv(AutoResolve):
	basedir: Path
	within_sample: Path = rel('sniffles_concat.self.tsv')
	across_patient: Path = rel('sniffles_concat.patient.tsv')
	def __post_init__(self): self._resolve()


@dataclass(frozen=True, slots=True)
class DataPaths(AutoResolve):
	# data dir & output dir (params when initialized)
	basedir: Path
	output: Path # This isn't really necessary here; it's part of cfg.
	# But, it's convenient to be able to access it from the DataPath object.
	
	# Sections
	ncbi: Ncbi = section(Ncbi)
	shc: Shc = section(Shc)
	isescan: IseScan = section(IseScan)
	lr_meta_sv: LrMetaSv = section(LrMetaSv)

	# Data files
	assembly_info: Path = rel('lr_assembly_info.csv') # Combined (NCBI nanopore flye_assembly_info + SH Plasmidsaurus computed assembly info)
	genomad: Path = rel('genomad.csv')
	contigs: Path = rel('contig_lengths.csv')
	gtdb: Path = rel('gtdbtk.csv')

	# Directories
	is_phylogeny: Path = rel('is_expansion_phylogeny') # dir within data containing relevant files
	timeline: Path = rel('shortread_is_timeline') # dir containing a subdirectory for each taxon, each of which contains a specific set of files
	
	def __post_init__(self): self._resolve()

