# Cop9_Signalosome Complex results

Directory for the main results of the modeling pipelines for the canonical and noncanonical CSN complex

## List of files and directories:

`results`   	 contains all results reported for both structures of the CSN complex
- IntegrativeStructure.CSN contains the subdirectories for the canonical structure of CSN complex (centroid and localization probability densities):
    * [Canonical CSN (DSSO+DHSO+BMSO)](./IntegrativeStructure.CSN/Structure_DSSO_DHSO_BMSO/) at a precision of 16Å.
    * [Subsample CSN (DHSO+BMSO)](./IntegrativeStructure.CSN/Structure_DHSO_BMSO/) at a precision of 22Å.
    * [Subsample CSN (DSSO+DHSO)](./IntegrativeStructure.CSN/Structure_DSSO_DHSO/) at a precision of 24Å.
    * [Subsample CSN (DSSO+BMSO)](./IntegrativeStructure.CSN/Structure_DSSO_BMSO/) at a precision of 27Å.
    * [Subsample CSN (DSSO)](./IntegrativeStructure.CSN/Structure_DSSO/) at a precision of 27Å.
    * [Subsample CSN (DHSO)](./IntegrativeStructure.CSN/Structure_DHSO/) at a precision of 29Å.
    * [Subsample CSN (BMSO)](./IntegrativeStructure.CSN/Structure_BMSO/) at a precision of 37Å.
- [Noncanonical CSN](./IntegrativeStructure.CSNn) contains the centroid and localization probability densities of the structure of the noncanonical CSN complex
- [Sampling Precision](./CSN_Sampling_Precision) contains the results of the exhaustiveness tests for the [canonical CSN complex](./Cross_Validations/DSSO_DHSO_BMSO	) and for each of the [subsamples](./CSN_Sampling_Precision/Cross_Validations/). The scripts are located in [Sampling_Precision_Scripts](./CSN_Sampling_Precision/Sampling_Precision_Scripts/). 

README is still under construction


_Author(s)_: Ilan E. Chemmama

_License_: [LGPL](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html).
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Publications_:
