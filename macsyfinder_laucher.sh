#!/bin/bash
for seq in Sequences/Annotations/*/*faa; do
    name=${seq##*/}
    for conj_type in typeF typeB typeC typeFATA typeFA typeG typeI typeT; do
        macsyfinder "$conj_type" \
                    -w 20 \
                    --db-type ordered_replicon \
                    -d Conjugation/Models/DEF \
                    -p Conjugation/Models/HMM \
                    --profile-suffix .hmm \
                    --sequence-db "$seq" \
                    -o Conjugation/${name%%.*}\_"$conj_type"
    done
done
