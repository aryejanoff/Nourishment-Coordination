The central file that determines the optimal rotation lengths (i.e., return periods between renourishments) for neighboring communities under coordination 
and non-coordination is "maincode.m". The first line in this script links an external script for parameter sensitivity analyses, and can be turned off if 
you wish to run the optimization for a specific set of parameters. Note that in order to do this, you must ensure that the parameters listed as inputs in 
the first line (i.e., "function [outputs]=maincode(inputs)") are turned on in the input parameters section of "maincode.m". This code requires access to the 
Parallel Computing Toolbox in Matlab.

Run file "basePV1_basePV2_sensitivity_analyses.m" to explore the model's sensitivity to different combinations of property values in each community. Sample 
output data from this model run is included in file "basePV1_basePV2_sensitivity_analyses_DATA.mat". Sample output data for the model run in figure 8b with
larger community sizes is also included in file "basePV1_basePV2_sensitivity_analyses_DATA_largerS.mat".

Run file "basePV1_basePV2_figures_behaviors_benefitofcoordination.m" to reproduce figures 5a-b, 5e, A1a-c in the paper (i.e., the emergent mode behaviors and 
the marginal importance of coordination for each property value combination in your regime space).

Run file "basePV1_basePV2_figures_efficiency.m" to reproduce figures 5c-d in the paper (i.e., physical nourishment efficiency for your regime space).

Run file "basePV1_basePV2_overundernourishment.m" to reproduce figures 5f-g, A1d-e in the paper (i.e., whether each community nourishes more/less than is 
economically optimal when choosing their strategy independently relative to their coordinated scheme).

Run file "model_rotation_ratio_analysis.m" to reproduce figures 8a-b in the paper (i.e., rotation length ratios for each paired community couplet). Envelopes
in this figure were created by drawing a box around the outer edges of the clustered model datapoints. Couple these results with the rotation length ratios 
from the field, contained in the spreadsheet with file path "Nourishment-Coordination/Field Data/NJ_Field_Data.xlsx".
