#!/bin/bash

blastn -db /genome_identification/hsv2/hsv2_irl_irs_id -query $1 -outfmt 6 > $1.irl_irs.txt
blastn -db /genome_identification/hsv2/hsv2_trl_id -query $1 -outfmt 6 > $1.trl.txt
blastn -db /genome_identification/hsv2/hsv2_trs_id -query $1 -outfmt 6 > $1.trs.txt
blastn -db /genome_identification/hsv2/hsv2_a_id -query $1 -outfmt 6 > $1.a.txt
blastn -db /genome_identification/hsv2/hsv2_us_id -query $1 -outfmt 6 > $1.us.txt
blastn -db /genome_identification/hsv2/hsv2_ul_id -query $1 -outfmt 6 > $1.ul.txt
