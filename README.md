# PAM - pellet ablation model

reference paper (https://iopscience.iop.org/article/10.1088/1741-4326/acb1c6/meta)

Changes from original python (OMFIT) version

1. do not use gEQDSK and most paramters from input.pam, take all of these from dd
2. make rho grid interpolation through IMASdd functions
3. all variables in MKS units 
4. different way to define a pellet position: first in x,y,z; then projection to r,z
