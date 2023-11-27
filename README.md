# Simultaneous convexification for the planar obnoxious facility location problem

This repo holds the global optima and source code for the algorithms described in the paper "Simultaneous convexification for the planar obnoxious facility location problem" submitted to the Journal of Global Optimization.

Instances are named `ofl_{N}_{M}{d}` where $N$ is the number of obnoxious facilities, $M$ is the number of communities, and $d \in$ `{i,ii}` refers to the minimum inter-facility squared distance $D$. If $d = $ `i`, then $D = 1/(2N)$, and if $d = $ `ii`, then $D = 1/N$. Global optima are provided as initial points in the GAMS files. The following termination criteria were used:
* `optca: 1e-5`
* `optcr: 1e-5`
* `reslim: 3600`

## Citation

Citation details will go here.