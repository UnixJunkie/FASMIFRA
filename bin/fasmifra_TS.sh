#!/bin/bash
#
# Copyright (C) 2024, Francois Berenger
# Tsuda laboratory, The University of Tokyo,
# 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
#
# FASMIFRA Thompson Sampling

#set -x # use -v on the CLI instead
set -u

# initial "population" of molecules from which fragments will be extracted
INIT_POP_FN="" # -p
# number of generated molecules at each iteration
PSIZE=""       # -s
# maximum number of iterations
MAXITER=""     # -i
# output directory; all created files will be in there
OUT=""         # -d
# scoring model file(s); one, or two coma-separated, filename(s)
MODEL=""       # -m

# CLI options handling
while getopts p:s:i:d:m:v opt; do
    case "${opt}" in
        p) INIT_POP_FN=${OPTARG}
           ;;
        s) PSIZE=${OPTARG}
           ;;
        i) MAXITER=${OPTARG}
           ;;
        d) OUT=${OPTARG}
           mkdir ${OUT} || exit 1
           ;;
        m) MODEL=${OPTARG}
           ;;
        v) set -x
           export OCAMLRUNPARAM='b'
           ;;
    esac
done
if [ "${INIT_POP_FN}" == "" ]; then echo "missing -p" && exit 1; fi
if [ "${PSIZE}" == "" ]; then echo "missing -s" && exit 1; fi
if [ "${MAXITER}" == "" ]; then echo "missing -i" && exit 1; fi
if [ "${OUT}" == "" ]; then echo "missing -d" && exit 1; fi
if [ "${MODEL}" == "" ]; then echo "missing -m" && exit 1; fi

###############################################################################
############ !!! EDIT this for your application !!! ###########################
###############################################################################

# WARNING: up to two (coma-separated) model files max
function score_molecules () {
    # one or two coma-separated files?
    if [[ "${MODEL}" =~ "," ]]; then
        echo 'predicting w/ two models in //...'
        M1=`echo ${MODEL} | cut -d',' -f1`
        M2=`echo ${MODEL} | cut -d',' -f2`
        # score in parallel
        molenc_gpr.py --load ${M1} -i $1 -o $2.1 2> $2.1.log &
        molenc_gpr.py --load ${M2} -i $1 -o $2.2 2> $2.2.log
        wait
        # average predictions to output file
        paste <(cut -f1,2 $2.1) <(cut -f1,2 $2.2) | \
            awk '($1 == $3){print $1"\t"($2+$4)/2.0}' > $2
    else
        # single model (standard case)
        molenc_gpr.py --load ${MODEL} -i $1 -o $2 2> $2.log
    fi
}

###############################################################################
###############################################################################
###############################################################################

cp ${INIT_POP_FN} ${OUT}/init_pop.smi

# activity statistics for the current molecular population
# (requires GNU datamash)
function stats () {
    cut -f2 $1 | datamash min 1 mean 1 sstdev 1 max 1 | \
        awk '{printf("\033[0;36m%.2f %.2f+/-%.2f %.2f\n\033[0m",$1,$2,$3,$4)}'
}

# check molecular diversity numbers
# if it is collapsing, then something is wrong
function mol_div_stats () {
    # compute cano SMILES
    molenc_smi2cansmi.py $1 > $1"_cano.smi"
    # monitor number of unique cano SMILES
    molenc_uniq -i $1"_cano.smi" -f 1 --force --in-RAM > $1"_cano_uniq.smi"
    # <parallel> --------------------------------------------------------------
    # monitor number of lead-like
    molenc_lead.py $1"_cano_uniq.smi" > $1"_cano_uniq_LL.smi" &
    # monitor number of drug-like
    molenc_drug.py $1"_cano_uniq.smi" > $1"_cano_uniq_DL.smi" &
    # monitor number of BM-scaffolds
    molenc_scaffold.py -i $1"_cano_uniq.smi" -o /dev/stdout 2>/dev/null | \
        cut -f3 | sort -u > $1"_cano_uniq_BM.smi"
    wait
    # </parallel> -------------------------------------------------------------
    printf "CANO/LEAD/DRUG/SCAF: %d %d %d %d\n" \
           `cat $1"_cano_uniq.smi"    | wc -l` \
           `cat $1"_cano_uniq_LL.smi" | wc -l` \
           `cat $1"_cano_uniq_DL.smi" | wc -l` \
           `cat $1"_cano_uniq_BM.smi" | wc -l`
}

# BRICS fragment it
~/src/FASMIFRA/bin/fasmifra_fragment.py --brics -i ${OUT}/init_pop.smi -o ${OUT}/init_pop_BRICS.smi 2> ${OUT}/BRICS.log

# generate 1st batch from it (no TS yet)
~/src/FASMIFRA/fasmifra -f -ufi -i ${OUT}/init_pop_BRICS.smi -n ${PSIZE} -o ${OUT}/gen_0.smi

# score 1st batch
score_molecules ${OUT}/gen_0.smi ${OUT}/gen_0.scores
echo "##### GEN 0" `stats ${OUT}/gen_0.scores`
mol_div_stats ${OUT}/gen_0.smi

# extract mu+/-s for TS
mean=`cut -f2 ${OUT}/gen_0.scores | datamash mean 1`
stdv=`cut -f2 ${OUT}/gen_0.scores | datamash sstdev 1`

# 1st TS iteration; -ig not possible
~/src/FASMIFRA/fasmifra -ufi -n ${PSIZE} -i ${OUT}/init_pop_BRICS.smi -o ${OUT}/gen_1.smi \
    --scores ${OUT}/gen_0.scores -mu ${mean} -s ${stdv} -og ${OUT}/gen_1.tsv
score_molecules ${OUT}/gen_1.smi ${OUT}/gen_1.scores
echo "##### GEN 1" `stats ${OUT}/gen_1.scores`
mol_div_stats ${OUT}/gen_1.smi

# subsequent TS iterations
# REMARK: we could use a patience=5 method watching the mean
#         as a smart stop condition...
for G in `seq 2 ${MAXITER}`; do
    GM1=`awk -v g=${G} 'BEGIN{print g-1}'`
    # generate population using TS
    ~/src/FASMIFRA/fasmifra -ufi -n ${PSIZE} -i ${OUT}/init_pop_BRICS.smi -o ${OUT}/gen_${G}.smi \
        --scores ${OUT}/gen_${GM1}.scores -ig ${OUT}/gen_${GM1}.tsv -og ${OUT}/gen_${G}.tsv
    # score it
    score_molecules ${OUT}/gen_${G}.smi ${OUT}/gen_${G}.scores
    echo "##### GEN" ${G} `stats ${OUT}/gen_${G}.scores`
    mol_div_stats ${OUT}/gen_${G}.smi
done
