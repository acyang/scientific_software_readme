###### tags: `software` `usage guide`

# NAMD quick usage guide

[![hackmd-github-sync-badge](https://hackmd.io/MvoYXS7ASB-DdNIvSmjPuw/badge)](https://hackmd.io/MvoYXS7ASB-DdNIvSmjPuw)

## Description
NAMD is a parallel molecular dynamics code designed for high-performance simulation of large biomolecular systems. NAMD uses the popular molecular graphics program VMD for simulation setup and trajectory analysis, but is also file-compatible with AMBER, CHARMM, and X-PLOR.

## Installation Status

| 功能           | Taiwania1  | Taiwania2     | Taiwania3  | TWCC container | TWCC VM      |
| :-------------: | :----------: | :-------------: | :----------: | :--------------: | :------------: |
| native cpu    | 2018.6<br>$$\surd$$  | 2.13<br>$$\star$$     | 2021.2<br>$$\star$$  | $$\triangle$$     | $$\triangle$$   |
| native gpu    | $$\times$$  | 2.13<br>$$\star$$     | $$\times$$  | $$\triangle$$     | $$\triangle$$   |
| CPU container | $$\triangle$$  | $$\triangle$$ | $$\triangle$$  | $$\triangle$$     | $$\triangle$$   |
| GPU container | 2.13<br>3.0<br>$$\star$$ | 2.13<br>3.0<br>$$\surd$$     | 2.13<br>3.0<br>$$\star$$ | $$\triangle$$     | $$\triangle$$   |

$$\surd \text{ : Tested} $$

$$\triangle \text{ : Not Ready} $$

$$\star \text{ : Untested} $$

$$\times \text{ : Not support} $$


## Installation Path
### Taiwania1
```/pkg/NAMD```
### Taiwania2
```/work/opt/ohpc/pkg/namd```
### Taiwania3
```/work/opt/ohpc/Taiwania3/pkg/chem/namd```

## Basic usage
### native app
```
source setnamd
namd2 +idlepoll +p 4 ${YOUR_INPUT}
```

### container
Use singularity to run application in container on HPC environment.
```
singularity run --nv -B ${PWD}:/host_pwd --pwd /host_pwd ${SIMG} namd2 ${YOUR_INPUT}
```

## Single Node Single GPU
 
* #### STEP1  Pull down the NAMD image
    ```docker pull nvcr.io/hpc/namd:3.0-alpha11```
* #### STEP2  Download the subm script
    ```export CWD=`pwd` ```
    ```cp /opt/ohpc/pkg/namd/submission/container/namd_snsg.sh $CWD```
* #### STEP3  Change your wallet account
    ```sed 's/JOB_ACCOUNT/XXXXXX/' namd_snsg.sh``` 
    
path 
* #### STEP4 Download simple example
    ```
    #!/bin/bash
    #SBATCH -J NAMD_Singlenode_SingleGPU_Job
    #SBATCH -A JOB_ACCOUNT
    #SBATCH --nodes=1
    #SBATCH --gres=gpu:1
    #SBATCH --ntasks-per-node=4
    #SBATCH --time=00:30:00
    #SBATCH --job-name=namd_job
    #SBATCH --output=namd_output.txt

    # module load
    module load singularity

    # set path 
    export APPROOT=$(pwd)
    INPUT="/home/yicheng0101/namd/apoa1/apoa1.namd"
    NAMD2="namd2 ${INPUT}"

    # set parameters
    export  NAMD_TAG=3.0-alpha11
    set -e; set -o pipefail
    GPU_COUNT=${1:-1}
    SIMG=${2:-"3.0-alpha11"}
    SINGULARITY="$(which singularity) exec --nv -B $(pwd):$APPROOT ${NAMD_TAG}.sif"
    echo $SINGULARITY

    # run the mpi job
    ${SINGULARITY} ${NAMD2}

    ```
    #### GPU 能用的數量
    ```
    GPU_COUNT=${1:-2}
    ```
    
    #### 指定要用哪幾張 GPUs
    ```
    export SINGULARITYENV_CUDA_VISIBLE_DEVICES=1,2
    ```

    XXXXXX 填入計畫帳號
* #### STEP4  Submit the job and execute
    ```sbatch namd_snsg.sh```
    
---

## Single Node Multiple GPU
* #### STEP1  Pull down the Gromacs image
    ```docker pull nvcr.io/hpc/namd:3.0-alpha11```
* #### STEP2  Download the subm script
    ```export CWD=`pwd` ```
    ```cp /opt/ohpc/pkg/namd/submission/container/namd_snmg.sh $CWD```
* #### STEP3  Change your wallet account
    ```sed 's/JOB_ACCOUNT/XXXXXX/' namd_snmg.sh``` 

* #### STEP4 Sets the number of processes to launch on each node
    ```
    ${SINGULARITY} ${NAMD_EXE} +ppn 2 +setcpuaffinity +idlepoll ${INPUT}
    ```
    ```
    #!/usr/bin/env bash
    #SBATCH -J GROMACS_Single_node_Multiple_GPUs_Job
    #SBATCH -A JOB_ACCOUNT
    #SBATCH --nodes=1
    #SBATCH --gres=gpu:1
    #SBATCH --ntasks-per-node=4
    #SBATCH --time=00:30:00
    #SBATCH --job-name=gromacs_job
    #SBATCH --output=gromacs_output.txt

    #module load
    module load singularity
    
    # Usage: ./singularity.sh <gpu count> <image name>
    set -e; set -o pipefail

    GPU_COUNT=${1:-2}
    SIMG=${2:-"3.0-alpha11"}

    echo "Downloading APOA1 Dataset..."
    INPUT="/home/yicheng0101/namd/apoa1/apoa1_nve_cuda.namd"
    export NAMD_EXE=namd2

    # Setting the GPU id
    TMP="0"
    for ((i=1;i<${GPU_COUNT};i++))
    do
      TMP="${TMP},${i}"
    done
    export   SINGULARITYENV_CUDA_VISIBLE_DEVICES=${TMP}

    # Run NAMD
    SINGULARITY="singularity exec --nv -B $(pwd):/home/yicheng0101/namd ${SIMG}.sif"
    NAMD2="namd2 ${INPUT}"

    echo "Running APOA1 example in ${SIMG} on ${GPU_COUNT} GPUS..."
    ${SINGULARITY} ${NAMD_EXE} +ppn 2 +setcpuaffinity +idlepoll ${INPUT}
    ```
Contributor: YI-CHENG HSIAO
    
