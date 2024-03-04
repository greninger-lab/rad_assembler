#!/bin/bash

blastn -db /genome_identification/hsv1/hsv1_irl_irs_id -query $1 -outfmt 6 > $1.irl_irs.txt
blastn -db /genome_identification/hsv1/hsv1_trl_id -query $1 -outfmt 6 > $1.trl.txt
blastn -db /genome_identification/hsv1/hsv1_trs_id -query $1 -outfmt 6 > $1.trs.txt
blastn -db /genome_identification/hsv1/hsv1_a_id -query $1 -outfmt 6 > $1.a.txt
blastn -db /genome_identification/hsv1/hsv1_us_id -query $1 -outfmt 6 > $1.us.txt
blastn -db /genome_identification/hsv1/hsv1_ul_id -query $1 -outfmt 6 > $1.ul.txt
