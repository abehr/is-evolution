
# safer bash defaults for all recipes
set shell := ["bash", "-eu", "-o", "pipefail", "-c"]
# auto-load environment variables from .env (so you don’t need to set PROJECT_CONFIG manually)
set dotenv-load := true

repo_root := source_directory() # dir of THIS file

# Set python to the venv-specific version, if/when it's installed. Otherwise use python3 in PATH
python := `if [ -x .venv/bin/python ]; then echo ".venv/bin/python"; else echo "python3"; fi`

# Check R
r := `if command -v Rscript >/dev/null; then echo "Rscript"; else echo "echo ⚠️ Rscript not found — R-based steps should be skipped"; fi`
