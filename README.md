<h2 align="center">Data analyses for project:" Brain-wide Resting-state fMRI Network Dynamics Elicited by Activation of Single Thalamic Input"
	
## Essential toolboxes

1. HMM-MAR-master: https://github.com/OHBA-analysis/HMM-MAR/
2. FASTR: https://fsl.fmrib.ox.ac.uk/eeglab/fmribplugin/
3. EEGLAB: https://sccn.ucsd.edu/eeglab/
4. spm12: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
5. DPABI： https://rfmri.org/DPABI/

## HMM models

We provide the HMM model trained on rsfMRI data from both OG– (n = 17, 90 trials) and OG+ (n = 11, 102 trials) conditions, saved as “Gamma_18_PCA70%.mat” and “hmm_18_PCA70%.mat”. This model was trained using PCA-reduced input, retaining 70% variance, and identified 18 hidden states.

The model was then applied to OG– and OG+ datasets separately, with the resulting state estimates saved as “Gamma_18_PCA70%_RS.mat” and “hmm_18_PCA70%_RS.mat” for OG–, and “Gamma_18_PCA70%_allOG.mat” and “hmm_18_PCA70%_allOG.mat” for OG+.

## HMM analysis codes

We provide the custom software codes used for HMM model training and estimation, computation and visualization of the state transition space, state decomposition, and generation of substate voxel-wise activation maps (without statistical thresholding). These procedures are implemented in the main script “Main_HMM_MAR.m”.
