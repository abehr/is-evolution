import os
import yaml
from pathlib import Path
from box import Box
from efm.data_paths import DataPaths

def resolve_nested_paths(d, root):
	for k, v in d.items():
		if k == 'project_root': continue # skip basedir
		if isinstance(v, dict):
			resolve_nested_paths(v, root)
		elif isinstance(v, str):
			# Convert to Path, and set relative to root if not an absolute path
			d[k] = Path(v).resolve() if v.startswith('/') else (root / v).resolve()


# ===== Find config file

# Find config file up to 3 layers back (prefer closer to current file vs. further back)
fp = Path(__file__).resolve() # current file path
cfg_names = ('config.yaml', 'config.yml') # valid names for config to search for
cfg_path_opts = [fp.parents[i] / fname for i in range(1,4) for fname in cfg_names]

# If env var exists, prefer that first (this can be any name).
env_cfg = os.environ.get('PROJECT_CONFIG')
if env_cfg: cfg_path_opts.insert(0, Path(env_cfg).expanduser())

for p in cfg_path_opts:
	if p.is_file(): break
if not p.is_file():
	raise FileNotFoundError('No config file found. Set $PROJECT_CONFIG or place config.yaml in a parent directory.')

# ===== Load config file

cfg = {}
with open(p) as f:
	cfg = Box(yaml.safe_load(f))


# ===== Check config directories

if not cfg.project_root.startswith('/'):
	raise ValueError(f'Error in config {p}: project_root {cfg.project_root} must be an absolute path.')

cfg.project_root = Path(cfg.project_root).resolve()

resolve_nested_paths(cfg, cfg.project_root)



# cfg.data_dir = Path(cfg.data_dir).resolve()
# cfg.output_dir = Path(cfg.output_dir).resolve()
# cfg.source = cfg.source # For now, don't process source any further. 

# Create output dir if it doesn't exist yet
cfg.output_dir.mkdir(parents=True, exist_ok=True)

for d in (cfg.project_root, cfg.data_dir, cfg.output_dir):
	if not d.is_dir():
		raise FileNotFoundError(f'Required directory {d} not found')


# ===== Load data

data = DataPaths(cfg.data_dir, cfg.output_dir)
