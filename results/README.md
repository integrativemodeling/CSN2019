# Cop9_Signalosome Complex results

Directory for the main results of the modeling pipelines for the canonical and noncanonical CSN complex

## List of files and directories:

`results`   	 contains all results reported for both structures of the CSN complex

- `CSN.cif` and `CSNn.cif` are the two CIFs file corresponding to the canonical and noncanonical structure of the CSN complex, respectively.

- `IntegrativeStructure.CSN` contains the subdirectories for the canonical structure of CSN complex (centroid and localization probability densities):
    * [Canonical CSN (DSSO+DHSO+BMSO)](./IntegrativeStructure.CSN/Structure_DSSO_DHSO_BMSO/) at a precision of 16Å.
    * [Subsample CSN (DHSO+BMSO)](./IntegrativeStructure.CSN/Structure_DHSO_BMSO/) at a precision of 22Å.
    * [Subsample CSN (DSSO+DHSO)](./IntegrativeStructure.CSN/Structure_DSSO_DHSO/) at a precision of 24Å.
    * [Subsample CSN (DSSO+BMSO)](./IntegrativeStructure.CSN/Structure_DSSO_BMSO/) at a precision of 27Å.
    * [Subsample CSN (DSSO)](./IntegrativeStructure.CSN/Structure_DSSO/) at a precision of 27Å.
    * [Subsample CSN (DHSO)](./IntegrativeStructure.CSN/Structure_DHSO/) at a precision of 29Å.
    * [Subsample CSN (BMSO)](./IntegrativeStructure.CSN/Structure_BMSO/) at a precision of 37Å.
    * [Sampling Precision](./sampling_precision_canonical) contains the results of the exhaustiveness tests for the [canonical CSN complex](./sampling_precision_canonical/DSSO_DHSO_BMSO) and for each of the [subsamples](./sampling_precision_canonical/sampling_precision_subsamples.tar). The plotting scripts are located in [sampling_precision_plotting_scripts](./sampling_precision_canonical/sampling_precision_plotting_scripts/). 

- `IntegrativeStructure.CSNn` contains the subdirectories for the noncanonical structure of CSN complex (centroid and localization probability densities):
    * [Noncanonical CSN](./IntegrativeStructure.CSNn/Structure_DSSO_DHSO_BMSO) at a precision of 22Å. 
    * [Sampling Precision](./IntegrativeStructure.CSNn/sampling_precision_noncanonical) contains the results of the exhaustiveness tests for the noncanonical CSN complex.
      
- `Mapping_XLs` contains the analysis of the cross-links on the different structures of the CSN complex.
    * [Canonical CSN](./Mapping_XLs/CSN_XL_Analysis) contains the results of the mapping of the canonical XLs to the canonical CSN complex.
    * [Noncanonical CSN](./Mapping_XLs/CSNn_XL_Analysis) contains the results of the mapping of the noncanonical XLs to the noncanonical CSN complex.
    * [Cross Mapping](./Mapping_XLs/CrossMapping) contains the mapping of all the cross-links to the crystal structure, canonical, and noncanonical structure of CSN. 

- `RMSD_ProteinLevel` contains the results of the structural differences among all the determined structure of CSN: canonical, subsampled canonical, noncanonical, crystal, bound to CRL1, and bound to CRL4. 
  
- `CSN_contact_frequency.zip` and `CSNn_contact_frequency.zip` are the files used to generate the contact frequency maps for the two CSN structures.
  
- `Chimera_sessions` contains different Chimera sessions to visualize the structures. 
  
  

README is still under construction


_Author(s)_: Ilan E. Chemmama

_License_: [LGPL](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html).
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Publications_:
