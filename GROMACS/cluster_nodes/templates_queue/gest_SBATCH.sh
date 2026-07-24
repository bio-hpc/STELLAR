#!/usr/bin/env bash
echo "#!/bin/sh" 										> $name_template_job
if [ "$renice" != "N/A" ];then
	echo "#SBATCH --nice=$renice"						>>$name_template_job 
fi
if [ "${GPU}" != "N/A" ];then
	echo "#SBATCH --gres=gpu:${GPU}" 					>>$name_template_job 
fi
echo "#SBATCH --output=${folder_out_jobs}${outJob}.out"	>>$name_template_job
echo "#SBATCH --error=${folder_out_jobs}${outJob}.err"	>>$name_template_job
echo "#SBATCH -J ${name_job}"							>>$name_template_job
if [ -n "${SLURM_PARTITION:-}" ]; then
	echo "#SBATCH --partition=${SLURM_PARTITION}"			>>$name_template_job
fi
echo "#SBATCH --time=${SLURM_TIME:-3-00:00:00}"			>>$name_template_job
echo "#SBATCH --cpus-per-task=${SLURM_CPUS_PER_TASK:-10}"	>>$name_template_job
echo "#SBATCH --nodes=1"								>>$name_template_job
echo "#SBATCH --ntasks=1"  							>>$name_template_job
if [ -n "${SLURM_QOS:-}" ]; then
	echo "#SBATCH --qos=${SLURM_QOS}"  					>>$name_template_job
fi
if [ -n "${SLURM_ACCOUNT:-}" ]; then
	echo "#SBATCH -A ${SLURM_ACCOUNT}"  					>>$name_template_job
fi

source ${path_cluster_nodes}templates_queue/codigoGR.sh
