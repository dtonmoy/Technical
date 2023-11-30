# -----
# Install a new kernelspec for python
# -----
python -m ipykernel install --user --name kernel_name

# -----
# Remove a kernelspec
# -----
jupyter kernelspec remove kernel_name

# -----
# Add a new kernelspec for R (inside R)
# -----
IRkernel::installspec(name = 'kernel_name', displayname = 'kernel_name')
