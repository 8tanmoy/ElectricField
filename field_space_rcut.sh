#!/bin/bash
for((ii=0;ii<75;ii++))
do
	mkdir rcut_$ii
	cp sub.sh rcut_$ii
	sed "s/_RCUTNO_/${ii}/g" field_space_rcut.py > rcut_${ii}/field_space_rcut.py
	cd rcut_$ii 
	qsub sub.sh
	cd ../
done