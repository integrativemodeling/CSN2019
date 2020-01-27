# Cop9_Signalosome Complex

Directory for the production scripts for the structure of the canonical and noncanonical CSN complex.

## List of files and directories:

`modeling_scripts`   	 contains all scripts used for modeling
- job.*.sub are the relevant submission scripts for a SGE cluster (Wynton at UCSF / QB3 Cluster at UCSF).
- smodeling.*.py are the modeling script for PMI. 

## Input for all the scripts

All inputs are located in the [data](../data) directory.

## Scripts used to generate the solution structures of the canonical and noncanonical CSN complex: 

- `job_cop9_dss_bms_dhs.sub` and `smodeling_dss_bms_dhs.py` are used in combination to generate the solution structure of the canonical CSN complex using BMSO, DHSO, and DSSO crosslinking data, excluded volution, and connectivity restraints.

- `job_cop9_dss_bms_dhs_plus.sub` and `smodeling_dss_bms_dhs_plus.py` are used in combination to generate the solution structure of the noncanonical CSN complex using BMSO, DHSO, and DSSO crosslinking data, excluded volution, and connectivity restraints.

## Scripts used to validate the structure of the canonical CSN complex by dataset-wide cross-validation

- `job_cop9_bms.sub` and `smodeling_bms.py` are used in combination to generate the structure of the canonical CSN complex with using BMSO crosslinking data, excluded volution, and connectivity restraints.  

- `job_cop9_dhs.sub` and `smodeling_dhs.py` are used in combination to generate the structure of the canonical CSN complex with using DHSO crosslinking data, excluded volution, and connectivity restraints.  

- `job_cop9_dss.sub` and `smodeling_dss.py` are used in combination to generate the structure of the canonical CSN complex with using DSSO crosslinking data, excluded volution, and connectivity restraints.  

- `job_cop9_dss_bms.sub` and `smodeling_dss_bms.py` are used in combination to generate the structure of the canonical CSN complex with using BMSO, DSSO crosslinking data, excluded volution, and connectivity restraints.

- `job_cop9_dss_dhs.sub` and `smodeling_dss_dhs.py` are used in combination to generate the structure of the canonical CSN complex with using DHSO, DSSO crosslinking data, excluded volution, and connectivity restraints.

- `job_cop9_bms_dhs.sub` and `smodeling_bms_dhs.py` are used in combination to generate the structure of the canonical CSN complex with using BMSO, DHSO crosslinking data, excluded volution, and connectivity restraints.


## Information

_Author(s)_: Ilan E. Chemmama

_License_: [LGPL](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html).
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Publications_:
