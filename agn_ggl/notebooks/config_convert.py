import os

file_base = 'GGL_AGN_UNIONS'

cmd = f'jupyter nbconvert --to script {file_base}.ipynb --stdout > {file_base}.py'

os.system(cmd)

# Run sciprt from command line via
# ipython script.py
