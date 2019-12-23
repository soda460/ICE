#!/bin/bash
i=1
for nc in "NC_000907" "NC_007146" "NC_009566" "NC_009567" "NC_014920" "NC_014922" "NC_016809" "NC_017452" "NC_017451" "NC_022356"; do
    printf -v j "%03d" $i
    echo HAIN"$j"
    prokka --kingdom Bacteria \
           --genus Haemophilus Sequences/Replicon/"$nc".fst \
           --outdir Sequences/Annotations/HAIN"$j" \
		   --proteins /data/ext4/dataDP/similarity_searches/db/card_original/protein_fasta_protein_homolog_model.fasta \
           --force \
           --locustag HAIN"$j" \
           --prefix HAIN"$j" \
           --cpus 0 >> LOG_prokka 2>&1
    i=$((i+1))
done
