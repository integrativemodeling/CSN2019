# Cop9_Signalosome Complex

Master data directory; used by all modeling scripts to generate the solutions of the CSN complex.


## List of files and directories:

- `data`   	 contains all relevant data
`CSN.fasta` contains all the sequences of the CSN complex (CSN1-CSN9).
`CSN.pdb` is the crystal structure of the canonical CSN complex (PDB code 4D10; [article](10.1038/nature13566))
`CSN1.pdb` is the complete model of the CSN1 subunit, obtained using MODELLER.
`CSN7B.pdb` is the structure of the CSN7B subunit isoform, modeled against subunit CSN7A isoform using MODELLER.

`BMS.Inter.csv` and `BMS.Intra.csv`: BMSO intermolecular and intramolecular cross-links for the canonical CSN complex
`DHS.Inter.csv` and `DHS.Intra.csv`: DHSO intermolecular and intramolecular cross-links for the canonical CSN complex
`DSS.Inter.csv` and `DSS.Intra.csv`: DSSO intermolecular and intramolecular cross-links for the canonical CSN complex

`BMS.Inter.P.csv` and `BMS.Intra.P.csv`: BMSO intermolecular and intramolecular cross-links for the CSN complex with CSN9
`DHSO.Inter.P.csv` and `DHSO.Intra.P.csv`: DHSO intermolecular and intramolecular cross-links for the CSN complex with CSN9
`DSS.Inter.P.csv` and `DSS.Intra.P.csv`: DSSO intermolecular and intramolecular cross-links for the CSN complex with CSN9

`topology_free.txt` and `topology_free_plus.txt` are the topology files for the canonical CSN complex and the CSN complex with CSN9, respectively.


## Information

_Author(s)_: Ilan E. Chemmama

_License_: [LGPL](http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html).
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Publications_:
