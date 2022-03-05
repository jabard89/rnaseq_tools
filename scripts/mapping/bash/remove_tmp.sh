#!/bin/bash
projdir="/home/jbard/beagle3-dadrummond/jbard/211216/STAR/output_220303"
for i in $(seq 30 99); do
	rm ${projdir}/JB${i}/*tmp
done
