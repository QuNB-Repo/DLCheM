# DLCheM on Compute Canada

**Login nodes** have access to the internet; **compute nodes** don't. That is why we adopt a two-stage process.

## Install the required dependencies in a virtual environment

The following **one-time** installation creates a virtual environment that contains all dependencies and can be used in compute jobs. 

1. Create a [virtual environment](https://docs.computecanada.ca/wiki/Python) in your home directory

```
mkdir -p ~/envs/dlchem
module load python/3.7
virtualenv --no-download ~/envs/dlchem/
```

2. Install the necessary libraries and [Jupyter scripts](https://docs.computecanada.ca/wiki/Jupyter) into that virtual environment

```
source ~/envs/dlchem/bin/activate
pip install --no-index torch ase tensorboardX h5py tqdm pytest
pip install schnetpack jupyter
echo -e '#!/bin/bash\nunset XDG_RUNTIME_DIR\njupyter notebook --ip $(hostname -f) --no-browser' > ~/envs/dlchem/bin/notebook.sh
chmod u+x ~/envs/dlchem/bin/notebook.sh
deactivate
```

## Running notebooks

After you have installed your virtual environment, you can use this environment to run notebooks. 

1. On a **login node**, navigate to the directory that contains your notebooks and load the virtual environment

```
cd ${directory containing notebooks}
module load python/3.7
source ~/envs/dlchem/bin/activate
```

2. Submit an [interactive job](https://docs.computecanada.ca/wiki/Running_jobs) with the required resources 

```
salloc --time=1:0:0 --ntasks=1 --cpus-per-task=2 --mem-per-cpu=1024M --account=def-stijn srun source ~/envs/dlchem/bin/notebook.sh
```

3. Following the `Connecting to Jupyter Notebook` instructions in the [Compute Canada Wiki](https://docs.computecanada.ca/wiki/Jupyter), create an SSH tunnel **locally**

```
sshuttle --dns -Nr <username>@<cluster>.computecanada.ca
```

and paste the required website address in your browser.
