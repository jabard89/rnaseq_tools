#!/bin/bash
projdir="/home/jbard/beagle3-dadrummond/jbard/211216/STAR/scripts"
for i in $(seq 30 99); do
	cp ${projdir}/220303_star_template.sh \
		${projdir}/sbatch/220303_star_JB${i}.sh; 
	sbatch ${projdir}/sbatch/map_JB${i}.sbatch
done
