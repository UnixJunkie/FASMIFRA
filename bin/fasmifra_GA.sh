#!/bin/bash
#
# Copyright (C) 2024, Francois Berenger
# Tsuda laboratory, The University of Tokyo,
# 5-1-5 Kashiwa-no-ha, Kashiwa-shi, Chiba-ken, 277-8561, Japan.
#
# FASMIFRA Genetic Algorithm

#set -x # use -v on the CLI instead
set -u

# initial "population" of molecules from which fragments will be extracted
INIT_POP_FN="" # -p
# number of generated molecules at each iteration
PSIZE=""       # -s
# PSIZE*GROWTH: size of the molecular population selected from
GROWTH=""      # -g
# maximum number of iterations
MAXITER=""     # -i
# output directory; all created files will be in there
OUT=""         # -d
# scoring model file(s); one, or two coma-separated, filename(s)
MODEL=""       # -m
VERBOSE=""     # -v

# CLI options handling
while getopts p:s:i:d:m:g:v opt; do
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
        g) GROWTH=${OPTARG}
           ;;
        v) set -x
           export OCAMLRUNPARAM='b'
           VERBOSE="TRUE"
           ;;
    esac
done
if [ "${INIT_POP_FN}" == "" ]; then echo "missing -p" && exit 1; fi
if [ "${PSIZE}" == "" ]; then echo "missing -s" && exit 1; fi
if [ "${MAXITER}" == "" ]; then echo "missing -i" && exit 1; fi
if [ "${OUT}" == "" ]; then echo "missing -d" && exit 1; fi
if [ "${MODEL}" == "" ]; then echo "missing -m" && exit 1; fi
if [ "${GROWTH}" == "" ]; then echo "missing -g" && exit 1; fi

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
    # LL <= 50 --> early stop
    if [ `cat $1"_cano_uniq_LL.smi" | wc -l` -le 50 ]; then
        echo "LL <= 50: early stop"
        exit 0
    fi
}

function show () {
    if [ "${VERBOSE}" != "" ]; then
        head -3 $1
    fi
}

# GA parameters: growth, elite proportion, diverse proportion
GROWN=`awk -v p=$PSIZE -v g=$GROWTH 'BEGIN{print int(p*g)}'`
DIVERSE=`awk -v p=$PSIZE 'BEGIN{print int(0.2*p)}'`
ELITE=`awk -v p=$PSIZE -v d=$DIVERSE 'BEGIN{print p-d}'`

cp ${INIT_POP_FN} ${OUT}/init_pop.smi

# BRICS fragment init. pop.
~/src/FASMIFRA/bin/fasmifra_fragment.py --brics -i ${OUT}/init_pop.smi -o ${OUT}/init_pop_BRICS.smi 2> ${OUT}/BRICS.log
show ${OUT}/init_pop_BRICS.smi

# generate 1st batch from it
~/src/FASMIFRA/fasmifra -f -pcb -i ${OUT}/init_pop_BRICS.smi -n ${GROWN} -o ${OUT}/curr_PCB.smi
show ${OUT}/curr_PCB.smi
~/src/FASMIFRA/bin/fasmifra_rm_cut_bonds.sh ${OUT}/curr_PCB.smi > ${OUT}/curr.smi
show ${OUT}/curr.smi

# score 1st batch
score_molecules ${OUT}/curr.smi ${OUT}/curr.scores
echo "##### GEN 0" `stats ${OUT}/curr.scores`
mol_div_stats ${OUT}/curr.smi

# generate next population
# diverse members
shuf -n $DIVERSE --random-source=/dev/urandom ${OUT}/curr_PCB.smi > ${OUT}/gen_1.smi
# elite members
paste ${OUT}/curr_PCB.smi ${OUT}/curr.scores > ${OUT}/curr.tsv
show ${OUT}/curr.tsv
# sort in place named (PCB) molecules and their scores
sort -n -r -k4 ${OUT}/curr.tsv -o ${OUT}/curr.tsv
show ${OUT}/curr.tsv
head -n $ELITE ${OUT}/curr.tsv | cut -f1,2 >> ${OUT}/gen_1.smi
show ${OUT}/gen_1.smi
~/src/FASMIFRA/bin/fasmifra_rm_cut_bonds.sh ${OUT}/gen_1.smi > ${OUT}/tmp.smi
show ${OUT}/tmp.smi
mol_div_stats ${OUT}/tmp.smi

# subsequent GA iterations
for G in `seq 2 ${MAXITER}`; do
    GM1=`awk -v g=${G} 'BEGIN{print g-1}'`
    # generate large population
    ~/src/FASMIFRA/fasmifra -f -pcb -i ${OUT}/gen_${GM1}.smi -n ${GROWN} -o ${OUT}/curr_PCB.smi
    show ${OUT}/curr_PCB.smi
    ~/src/FASMIFRA/bin/fasmifra_rm_cut_bonds.sh ${OUT}/curr_PCB.smi > ${OUT}/curr.smi
    show ${OUT}/curr.smi
    # score it
    score_molecules ${OUT}/curr.smi ${OUT}/curr.scores
    show ${OUT}/curr.scores
    echo "##### GEN "$G `stats ${OUT}/curr.scores`
    # generate next population
    # diverse members
    shuf -n $DIVERSE --random-source=/dev/urandom ${OUT}/curr_PCB.smi > ${OUT}/gen_${G}.smi
    # elite members
    paste ${OUT}/curr_PCB.smi ${OUT}/curr.scores > ${OUT}/curr.tsv
    show ${OUT}/curr.tsv
    # sort in place named (PCB) molecules and their scores
    sort -n -r -k4 ${OUT}/curr.tsv -o ${OUT}/curr.tsv
    show ${OUT}/curr.tsv
    head -n $ELITE ${OUT}/curr.tsv | cut -f1,2 >> ${OUT}/gen_${G}.smi
    show ${OUT}/gen_${G}.smi
    ~/src/FASMIFRA/bin/fasmifra_rm_cut_bonds.sh ${OUT}/gen_${G}.smi > ${OUT}/tmp.smi
    show ${OUT}/tmp.smi
    mol_div_stats ${OUT}/tmp.smi
done
