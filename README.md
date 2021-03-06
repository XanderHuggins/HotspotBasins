## Code for: Hotspots for social and ecological impacts from freshwater stress and storage loss <br/>
[![DOI](https://zenodo.org/badge/419507061.svg)](https://zenodo.org/badge/latestdoi/419507061)

This is the code repository associated with the paper Huggins et al. Hotspots for social and ecological impacts from freshwater stress and storage loss. *Nature Communications* (2022). https://doi.org/10.1038/s41467-022-28029-w. <br/>

**Repository structure:** <br/>
*Note: Individual scripts are annotated and a description is provided at the top of each script.*  <br/>
The order of the folders below is the required workflow to reproduce any results from the study. Code numberings within folders, e.g. `1-...`, `2-...`, etc., indicate either the necessary or suggested order to sequence the codes. 

* **R/** -- 
    * **setup/** -- scripts that import necessary packages and set common plotting arguments
    * **udf/** -- scripts with custom functions and are described at the top of each script
    * **preprocessing/** -- scripts that spatially harmonize and otherwise prepare data for analysis 
    * **analysis/** -- scripts that perform core analyses of this study, including uncertainty and sensitivity analyses
    * **plotting+stats/** -- scripts that create all plots and calculate summary statistics reported in the study
