#!/bin/bash
for((ii=0;ii<=1000;ii++))	#1000
do
	mkdir t_$ii
	cp sub.sh t_$ii
	sed "s/_TNO_/${ii}/g" field_space_projection.py > t_${ii}/field_space_projection.py
	cd t_$ii 
	qsub sub.sh
	cd ../
done