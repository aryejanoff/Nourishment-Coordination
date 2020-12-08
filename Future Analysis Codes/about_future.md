The central file that determines the optimal rotation lengths (i.e., return periods between renourishments) for neighboring communities under coordination 
and non-coordination is "maincode.m". The first line in this script links an external script for parameter sensitivity analyses, and can be turned off if 
you wish to run the optimization for a specific set of parameters. Note that in order to do this, you must ensure that the parameters listed as inputs in 
the first line (i.e., "function [outputs]=maincode(inputs)") are turned on in the input parameters section of "maincode.m". This code requires access to the 
Parallel Computing Toolbox in Matlab.

Run file "future_erosion_sandcost_sensitivity_analyses.m" to explore the model's sensitivity to different combinations of background erosion rates and sand
costs in each community. Sample output data from this model run is included in file "future_erosion_sandcost_DATA.mat".

Run file "future_erosion_sandcost_figures_behaviors.m" to reproduce figures 10a-b, A2a-c in the paper (i.e., the emergent mode behaviors and 
the marginal importance of coordination for each erosion rate/sand cost combination in your regime space).

Run file "future_erosion_sandcost_figures_efficiency.m" to reproduce figures 10d-e in the paper (i.e., physical nourishment efficiency for your regime space).
