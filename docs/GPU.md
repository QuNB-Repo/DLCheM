<p align="center">
<img src="../media/hpc-gpu.png" width="400">
</p>

This file is intended to document the idiosyncrasies of the Compute Canada HPC systems that provide access to GPUs. We assume that you know [how to run jobs](https://docs.computecanada.ca/wiki/Running_jobs) that use [Python](https://docs.computecanada.ca/wiki/Python) and (general) [GPUs](https://docs.computecanada.ca/wiki/Using_GPUs_with_Slurm). Note that [specific performance issues](https://docs.computecanada.ca/wiki/AI_and_Machine_Learning) are important for Machine Learning, which are also introduced in the respective [Compute Canada tutorials](https://docs.computecanada.ca/wiki/Tutoriel_Apprentissage_machine/en).

# General Compute Canada infrastructure

The Compute Canada GPU infrastructure is documented [here](https://docs.computecanada.ca/wiki/Using_GPUs_with_Slurm#Available_hardware), alongside the [Slurm type specifiers](https://docs.computecanada.ca/wiki/Running_jobs). According to [this website](https://www.microway.com/knowledge-center-articles/comparison-of-nvidia-geforce-gpus-and-nvidia-tesla-gpus/), useful performance metrics are:

| Type of GPU | TFLOPS64<sup>[1](#myfootnote1)</sup> | TFLOPS16<sup>[2](#myfootnote2)</sup> | TensorFlops<sup>[3](#myfootnote3)</sup> | Mem Quantity (GB) | Mem performance (GB/s) | Comments |
| --- | --- | --- | --- | --- | --- | --- |
| P100 | 4.7 | 18.7 | N/A | 12/16<sup>[4](#myfootnote4)</sup> | 549/732 |
| V100 | **7** | 28 | 112 | 16/32<sup>[4](#myfootnote4)</sup> | 900 |
| T4 | 0.25 | 16.2 | 65 | 16<sup>[4](#myfootnote4)</sup> | 320 |
| K80 | 1.87 | N/A | N/A | 24 | 480 | 
| RTX6000 | 0.5 | **32.6** | **130.5** | **24**<sup>[4](#myfootnote4)</sup> | 624 | Currently only available on Siku |

Bottom line: use `V100`s or `RTX6000` when possible.

<a name="myfootnote1">1</a>: FP64 64-bit (Double Precision) Floating Point Calculations.

<a name="myfootnote2">2</a>: FP16 16-bit (Half Precision) Floating Point Calculations.

<a name="myfootnote3">3</a>: combines a multiply of two FP16 units (into a full precision product) with a FP32 accumulate operation.

<a name="myfootnote4">4</a>: Tesla/Quadro Unified Memory allows GPUs to share each other’s memory.

# Siku 

<p align="center">
<img src="../media/siku.jpeg" width="200">
</p>

## Getting SSH access

You will need to **manuall** install your public SSH key on the login nodes by copying the contents of your ~/.ssh/*.pub to the file ~/.ssh/authorized_keys on the Siku login node. If needed, create the `authorized_keys` file yourself 

```bash
touch ~/.ssh/authorized_keys
```

## Check that you can run on the provided GPUs

Submit an interactive job


```bash
salloc --time=1:00:00 --cpus-per-task=10 --mem=46GB --gres=gpu:rtx6000:1 
--account=ctb-stijn(or def-stijn) --partition=all_gpus
```
note that if you choose ctb-stijn account you are limited to 80 cores 
and 6 gpus but have much higher priority and can run for up to 72 hours

if you choose def-stijn you are no longer restricted to run only 80 
cores and 6 gpus, but you will have lower priority and can run up to 24 
hours

load the needed modules when your job has started

```
module load python cuda cudnn
```

and create a [Python virtual environment](https://docs.computecanada.ca/wiki/Python)

```
cd $SLURM_TMPDIR
virtualenv --no-download env
source $SLURM_TMPDIR/env/bin/activate
pip install --no-index --upgrade pip
pip install --no-index -r requirements.txt
```

where `requirements.txt` can be found [here](requirements.txt). After your virtual environment has been created, run the following commands inside `ipython`.

```python
import torch
torch.cuda.is_available()
torch.cuda.get_device_name()
```

The last command should return `'Quadro RTX 6000'`.

Again, as before, to make jupyter notebook executible, run:

```bash
echo -e '#!/bin/bash\nunset XDG_RUNTIME_DIR\njupyter notebook --ip 
$(hostname -f) --no-browser' > ~/slurm/env/bin/notebook.sh
```
now you are ready to run the executable

```
chmod u+x ~/slurm/env/bin/notebook.sh
```

now you can run using:

```
srun ~/slurm/env/bin/notebook.sh
```
