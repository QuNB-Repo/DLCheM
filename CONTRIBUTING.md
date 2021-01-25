# DLCheM on Compute Canada

**Login nodes** have access to the internet; **compute nodes** don't. That is why we adopt a two-stage process.

## Admin: installation

Create a [virtual environment](https://docs.computecanada.ca/wiki/Python) in our [shared envs](https://docs.computecanada.ca/wiki/Sharing_data)

```
mkdir ~/projects/def-stijn/envs/dlchem/
module load python/3.7
virtualenv --no-download ~/projects/def-stijn/envs/dlchem/
```

Install the necessary libraries into that virtual environment

```
source ~/projects/def-stijn/envs/dlchem/bin/activate
pip install --no-index torch ase tensorboardX h5py tqdm pytest
pip install schnetpack jupyter
echo -e '#!/bin/bash\nunset XDG_RUNTIME_DIR\njupyter notebook --ip $(hostname -f) --no-browser' > ~/projects/def-stijn/envs/dlchem/bin/notebook.sh
chmod ug+x ~/projects/def-stijn/envs/dlchem/bin/notebook.sh
deactivate
```

## Users: running notebooks

Load the Python module and the virtual environment on a **login node**

```
module load python/3.7
source ~/projects/def-stijn/envs/dlchem/bin/activate
```

Then submit an interactive job with the required resources
```
salloc --time=1:0:0 --ntasks=1 --cpus-per-task=2 --mem-per-cpu=1024M --account=def-stijn srun ~/projects/def-stijn/envs/dlchem/bin/notebook.sh
```

Following the `Connecting to Jupyter Notebook` instructions in the [Compute Canada Wiki](https://docs.computecanada.ca/wiki/Jupyter), create an SSH tunnel **locally**
```
sshuttle --dns -Nr <username>@<cluster>.computecanada.ca
```

and paste the required website address in your browser.
