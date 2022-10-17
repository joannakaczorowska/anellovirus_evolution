#!/bin/bash

##This is a script for MEME
hyphy meme --alignment input_hyphy_sep2022/consensus_S2_1_FEL.fasta \
--tree input_hyphy_sep2022/consensus_S2_1_FEL_MAFFT_tree.nwk \
--output input_hyphy_sep2022/consensus_S2_1_MEME.json
