# -----
# Create a new conda environment inside a directory with prefix
# -----
conda create --prefix ./env pandas numpy matplotlib scikit-learn

# -----
# Create a new conda environment with name
# -----
conda create --name env_name pandas numpy matplotlib scikit-learn

# -----
# Remove a conda environment
# -----
conda env remove -–prefix "conda_env"

# -----
# clone an existing conda repository using prefix
# -----
mamba create --prefix new_env --clone path/to/file/env1

# -----
# Remove a conda environment based on conda repository
# -----
conda remove --name "conda_env" --all

