#!/bin/bash

blastn -db /genome_identification/cmv/cmv_irl_irs_id -query $1 -outfmt 6 > $1.irl_irs.txt
blastn -db /genome_identification/cmv/cmv_hyper_id -query $1 -outfmt 6 > $1.hyper.txt
blastn -db /genome_identification/cmv/cmv_trl_id -query $1 -outfmt 6 > $1.trl.txt
blastn -db /genome_identification/cmv/cmv_trs_id -query $1 -outfmt 6 > $1.trs.txt
blastn -db /genome_identification/cmv/cmv_trs1_id -query $1 -outfmt 6 > $1.trs1.txt
blastn -db /genome_identification/cmv/cmv_pretrs1_id -query $1 -outfmt 6 > $1.pretrs1.txt
blastn -db /genome_identification/cmv/cmv_a_id -query $1 -outfmt 6 > $1.a.txt
blastn -db /genome_identification/cmv/cmv_us_id -query $1 -outfmt 6 > $1.us.txt
blastn -db /genome_identification/cmv/cmv_ussub_id -query $1 -outfmt 6 > $1.ussub.txt
