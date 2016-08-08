#!/bin/sh
# quick utility to detect unused variables in the header files
# it is far from perfect, but for a quick and easy tool it does its job

for header in `ls *.h`
do
	for T in "double" "int" "bool" "vec3d" "Averager" "StressTensor"
	do
		for v in `grep $T $header | grep ";" | awk '{print $2}' | grep ";" | sed -e "s/;//g" | grep -v "(" | sed -e "s/\*//g"` # quick and dirty parsing of members of type T in header
		do
			if [ `grep "$v" *.cpp *.h | wc -l` -lt "2" ] # only one match: in the header
			then
				echo $header ": " $T $v
			fi
		done
	done
done
