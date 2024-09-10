
#!/bin/bash           # the shell language when run outside of the job scheduler
#                     # lines starting with #$ is an instruction to the job scheduler
#$ -S /bin/bash       # the shell language when run via the job scheduler [IMPORTANT]
#$ -cwd               # job should run in the current working directory
#$ -N sim_example_4W52 # set the name of the job
#$ -j y               # STDERR and STDOUT should be joined
#$ -o experiments/example_4W52/simulations/       # set the destination for the STDOUT file
#$ -l mem_free=4G     # job requires up to 1 GiB of RAM per slot
#$ -l scratch=2G      # job requires up to 2 GiB of local /scratch space
#$ -l h_rt=2:00:00    # job requires up to this many hours of runtime
#$ -r y               # if job crashes, it should be restarted
#$ -q gpu.q           # use the gpu queue

## load the required modules
module load CBI
module load miniconda3/4.12.0-py39
module load Sali
module load cuda/10.0.130

## print start time:
date

## print all cuda devices:
echo $CUDA_VISIBLE_DEVICES
## make sure we only use the assigned GPU:
export CUDA_VISIBLE_DEVICES=$SGE_GPU
echo $CUDA_VISIBLE_DEVICES

EXPERIMENT_DIR=$1

## Run the simulation
apptainer exec --nv $EASYMD_CONTAINER bash -c \
". /opt/conda/etc/profile.d/conda.sh && \
conda activate /opt/conda/envs/easyMD && \
python /easyMD/autorun.py $EXPERIMENT_DIR"

## End-of-job summary, if running as a job
[[ -n "$JOB_ID" ]] && qstat -j "$JOB_ID"  # This is useful for debugging and usage purposes,
                                        # e.g. "did my job exceed its memory request?"

## print end time:
date
    