# RAGA Dataset

# Introduction

This repository contains the RAGA dataset. It is an affective a labelled dataset which use racing games (in VR and with a standard desktop) as stimuli. In particular, it involves physiological data (5 facial EMG, ECG, EDA, and Respiration) alligned with the ground truth. The latter has been expressed in core affect, i.e., valence and arousal data. It has been labeled through a self-assessment over all the game session. For more details see [An Empirical Study of Players’ Emotions in VR Racing Games Based on a Dataset of Physiological Data](https://link.springer.com/article/10.1007/s11042-019-08585-y) by Granato et al.

# Citation

Granato, M., Gadia, D., Maggiorini, D. et al. An empirical study of players’ emotions in VR racing games based on a dataset of physiological data. Multimed Tools Appl (2020). https://doi.org/10.1007/s11042-019-08585-y

Bibtex:

```tex
@article{granato2020empirical,
  title={An empirical study of players’ emotions in VR racing games based on a dataset of physiological data}, 
  DOI={10.1007/s11042-019-08585-y}, 
  journal={Multimedia Tools and Applications}, 
  author={Granato, Marco and Gadia, Davide and Maggiorini, Dario and Ripamonti, Laura A.}, 
  year={2020}, 
  month={Feb}
}
```
  

# Directory Structure

The directories available in this repository follows this structure:

- "./" The source code of matlab file to analyse the data
- "./Functions" contains the functions for the matlab source code
- "./Data/Examples" The test files
- "./Data/Real" the real data 
- "./MatlabDataset" where are stored the .mat files
- "./Results/" where are stored the files of the results


NB: The data have been stored with the following directory structure: id/vr; id/novr (eg 01/vr;01/novr);



