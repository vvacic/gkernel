#!/bin/bash

DATASET=LYS_3LCK

PDB_DIR=./data
GRAPH_DIR=./data
OUT_DIR=./data

CONNECTION=C_ALPHA
RADIUS=6
SPHERE=30
REDUCTION=BLOSUM_15

date

echo
echo Generating graph files for each residue of interest...

for RESIDUE in 231 246 269 273 276 293 329 335 340 379 401 405 420 478; do
echo $RESIDUE
./parse_pdb -c A -m ${CONNECTION} -d ${RADIUS} -s ${SPHERE} -C -r ${RESIDUE} -o ${OUT_DIR}/3LCK_A_${RESIDUE}.graph ${PDB_DIR}/3LCK.pdb
done

echo
echo Generating kernel matrix...

./make_kernel \
-p ${GRAPH_DIR}/${DATASET}.positives \
-n ${GRAPH_DIR}/${DATASET}.negatives \
-a $GRAPH_DIR \
-r $REDUCTION \
-N \
-V \
-s ${OUT_DIR}/${DATASET}_${REDUCTION}.svml

./make_kernel \
-p ${GRAPH_DIR}/${DATASET}.positives \
-n ${GRAPH_DIR}/${DATASET}.negatives \
-a $GRAPH_DIR \
-r $REDUCTION \
-N \
-V \
-k ${OUT_DIR}/${DATASET}_${REDUCTION}.dat

echo

date

