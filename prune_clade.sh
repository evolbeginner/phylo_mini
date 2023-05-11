#! /bin/bash


######################################################################
treefile=''
outgrp=''
is_v=false


######################################################################
while [ $# -gt 0 ]; do
	case $1 in
		-t)
			treefile=$2
			shift
			;;
		--outgroup|--outgrp)
			outgrp=`echo $2 | sed 's/,/ /g'`
			shift
			;;
		-v)
			is_v=true
			;;
		*)
			echo "wrong argu $1" >&2
			exit 1
	esac
	shift
done


######################################################################
outgrp_species=`nw_clade $treefile $outgrp |nw_distance -n - | cut -f1`


if [ $is_v == true ]; then
	nw_prune -v $treefile $outgrp_species
else
	nw_prune $treefile $outgrp_species
fi


