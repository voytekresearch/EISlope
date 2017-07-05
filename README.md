#Inferring synaptic excitation/inhibition balance from field potentials
Analysis and plotting code from Gao et al., 2017 NeuroImage
(http://www.sciencedirect.com/science/article/pii/S1053811917305621)

This repo contains scripts (./scripts) used to create figures shown in the above referenced publication. All non-MATLAB dependencies are included in ./functions. 
MATLAB version: R2014b

Figure 1: ./scripts/gao2017_PSDsim.m
Figure 2: ./scripts/gao2017_ca1depth.m & EI_models.ipynb
Figure 3: ./scripts/gao2017_ca1phase.m
Figure 4: ./scripts/gao2017_tycho.m

The simulation figures (Figure 1) can be reproduced exactly by running the script. For figures 2-4, some example data are included in example_data, where a single session of rat CA1 LFP data and one channel of monkey ECoG are included. For obvious reasons (space), the full dataset is not uploaded on Github. 

One can find these publicly available datasets at:
rat: https://crcns.org/data-sets/hc/hc-2/about-hc-2
monkey (propofol): http://neurotycho.org/anesthesia-and-sleep-task