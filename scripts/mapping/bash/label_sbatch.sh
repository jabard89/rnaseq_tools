#!/bin/bash
for i in $(seq 30 99);	do
	sed s/'${BASHINDEX}'/${i}/g map_template.sbatch > sbatch/map_JB${i}.sbatch
done
