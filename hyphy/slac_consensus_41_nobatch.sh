#!/bin/bash

##This is a script for SLAC
hyphy slac --alignment input_hyphy_sep2022/consensus_41_FEL.fasta \
--tree input_hyphy_sep2022/consensus_41_FEL_MAFFT_tree.nwk \
--output input_hyphy_sep2022/consensus_41_SLAC.json
