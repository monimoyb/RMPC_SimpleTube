# RMPC_SimpleTube
A simple Robust MPC for Linear Systems with Parameric and Additive Uncertainty. These codes produce the results of our proposed algorithm, shown in the paper https://arxiv.org/abs/2103.12351. The Region of Attraction (ROA) of our controller (yellow set) contains about 98% of the ROA of Tube MPC (grey set, https://www.sciencedirect.com/science/article/pii/S0005109803002838?via%3Dihub, Section 5). The yellow ROA is about 4% larger in volume.

<img width="889" alt="roa" src="https://user-images.githubusercontent.com/12418616/113952864-685cce80-97cb-11eb-8601-2cec242980db.png">

The following table of computation times demonstrates that our approach is up to 15x faster at online control computations compared to the aforementioned Tube MPC.
![time](https://user-images.githubusercontent.com/12418616/113953069-d608fa80-97cb-11eb-87dc-bd13e1152b14.png)

Additionally, even with an open-loop policy sequence (i.e., without re-solving the MPC problem), we obtain a safe set (yellow set), which is about 12x larger in volume than the ROA of the System Level Synthesis based Constrained LQR (blue set, https://arxiv.org/pdf/1809.10121.pdf, Section 2.3).
<img width="911" alt="roa_sls" src="https://user-images.githubusercontent.com/12418616/113953303-62b3b880-97cc-11eb-9fc3-57bc0d422495.png">
