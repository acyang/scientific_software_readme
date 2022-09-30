# GROMACS quick usage guide
## Description
GROMACS is a versatile package to perform molecular dynamics, i.e. simulate the Newtonian equations of motion for systems with hundreds to millions of particles and is a community-driven project. Contributions are welcome in many forms, including improvements to documentation, patches to fix bugs, advice on the forums, bug reports that let us reproduce the issue, and new functionality. 
## Single Node Single GPU
cp 
* #### STEP1  Pull down the Gromacs image
    ```docker pull nvcr.io/hpc/gromacs:2022.1```
* #### STEP2  Download the subm script
    ```export CWD=`pwd` ```
    ```cp /opt/ohpc/pkg/gromacs/container/test/gromacs_single_node_single_gpu.sh $CWD```
* #### STEP3  Change your wallet account
    ```sed 's/JOB_ACCOUNT/XXXXXX/' gromacs_snsg.sh``` 
        
    XXXXXX 填入計畫帳號
* #### STEP4  Submit the job and execute
    ```sbatch gromacs_snsg.sh```

dataset name
* #### STEP5 Prepare your input files - Run a simple testcase
    * ##### STEP1 Download the testcase
        ```wget -c ftp://ftp.gromacs.org/pub/benchmarks/${DATA_SET}.tar.gz```
    * ##### STEP2 Decompression the file 
        ```tar -xvf water_GMX50_bare.tar.gz ```
    * ##### STEP3 Choose one testcase
        ```cd water-cut1.0_GMX50_bare/1536/```
    * ##### STEP4 Change the format
        ```singularity gmx grompp -f pme.mdp``` 
    * ##### STEP5 Change the testcase - Run benchmark
        ```
        singularity gmx mdrun \
                     -ntmpi ${GPU_COUNT} \
                     -nb gpu \
                     -ntomp ${OMP_NUM_THREADS} \
                     -pin on \
                     -v \
                     -noconfout \
                     -nsteps 5000 \
                     -s topol.tpr
        ```
    * ```singularity gmx mdrun```
    * ```-ntmpi ${GPU_COUNT}```
    * ```-nb gpu```
    * ```-ntomp ${OMP_NUM_THREADS}```
    * ```-pin on```
    * ```-noconfout```
    * ```-nsteps 5000 ```
    * ```-s topol.tpr ``` Test case 產生的tpr 檔

    #### gromacs_snsg.sh
   ```
    #!/bin/bash

    #SBATCH -J GROMACS_Singlenode_SingleGPU_Job
    #SBATCH -A JOB_ACCOUNT
    #SBATCH --nodes=1
    #SBATCH --gres=gpu:1
    #SBATCH --ntasks-per-node=4
    #SBATCH --time=00:30:00
    #SBATCH --job-name=gromacs_job
    #SBATCH --output=gromacs_output.txt

    #module load
    module load singularity

    # Initialize the path
    export APP_ROOT=$PWD
    cd $APP_ROOT

    # Script arguments
    GPU_COUNT=${1:-1}
    SIMG=${2:-"${PWD}/gromacs-2022_1.sif"}

    # Set number of OpenMP threads
    export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}

    # Singularity will mount the host PWD to /host_pwd in the container
    SINGULARITY="singularity run --nv -B ${PWD}:/host_pwd --pwd /host_pwd ${SIMG}"

    # Prepare benchmark data
    ${SINGULARITY} gmx grompp -f pme.mdp

    # Run benchmark
    ${SINGULARITY} gmx mdrun \
                     -ntmpi ${GPU_COUNT} \
                     -nb gpu \
                     -ntomp ${OMP_NUM_THREADS} \
                     -pin on \
                     -v \
                     -noconfout \
                     -nsteps 5000 \
                     -s topol.tpr


Contributor: YI-CHENG HSIAO