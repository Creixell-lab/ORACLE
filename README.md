# ORACLE: ORigins and Context-independent Analysis of CDMS Libraries for Epistasis

This is the official repository of ORACLE for the paper:
Identifying Protein Superbinders and Molecular Determinants of Epistasis with Combinatorial Deep Mutational Scanning (CDMS) Libraries

Mingxuan Jiang<sup>1</sup>, Mohan Sun<sup>1</sup>, Nuo Cheng<sup>1</sup>, Mihkel Örd<sup>1</sup>, Teresa L. Augustin<sup>1</sup>,  
Allyson Li<sup>2</sup>, Neel H. Shah<sup>2</sup>, Jesse Rinehart<sup>3,4</sup>, Helen R. Mott<sup>5</sup>, Pau Creixell<sup>1</sup><sup>*</sup>

## Affiliations

<sup>1</sup> Cancer Research UK Cambridge Institute, University of Cambridge, Li Ka Shing Centre, Robinson Way, Cambridge CB2 0RE, UK  
<sup>2</sup> Department of Chemistry, Columbia University, 3000 Broadway, New York, NY 10027, USA  
<sup>3</sup> Department of Cellular & Molecular Physiology, Yale School of Medicine, New Haven, CT 06520, USA  
<sup>4</sup> Systems Biology Institute, Yale University, West Haven, CT 06516, USA  
<sup>5</sup> Department of Biochemistry, University of Cambridge, Tennis Court Road, Cambridge CB2 1QW, UK  

<sup>*</sup> Corresponding author

Code for WT independent Epistasis Calculations
![image](https://github.com/user-attachments/assets/828ab39c-3dbe-401e-9115-7d4a017cb8cb)





[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Creixell-lab/epystasis/blob/main/Epistasis_20250811.ipynb)

Google colab file is at:
https://colab.research.google.com/github/Creixell-lab/epystasis/blob/main/Epistasis_20250811.ipynb
This colab file allows an interactive parsing of the underlying logic and pipeline of the approach.


ORACLE is a WT-agnostic, Walsh-Hadamard based framework that can treats phenotype-genotype datasets using a CDMS (combinatorial deep mutational scanning) approach with interpretable epistasis modelling to generate 1st order, 2nd order, 3rd order and higher order of epistasis coefficients. With this coefficients, we can quantify protein interaction landscapes and show interactions that can be explained with single mutations (1st order) and interactions that rely on cumulative, additive interactions of amino acids. Using these epistasis coefficients can also help to predict out of sample phenotypes for missing or unsampled variants.

## Installation guide

ORACLE epistasis can also be run locally using the following pip installation guide
...
