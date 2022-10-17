#!/bin/bash

##This is a script for slac
hyphy slac --alignment input_hyphy_sep2022/consensus_23_FEL.fasta \
--tree input_hyphy_sep2022/consensus_23_FEL_MAFFT_tree.nwk \
--output input_hyphy_sep2022/consensus_23_SLAC.json
