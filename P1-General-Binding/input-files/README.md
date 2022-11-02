Sample files required to run a MELD run (from 1A74 system):
1. restraints.py
Restraint functions are defined in this file.
2. setupMELD.py
All inputs including simulations details, restraints and other options are specified here.
3. protein-contacts.dat
Tertiary structure contact list needed to prevent denaturation.
4. sequence.dat
Sequence of DNA needed to generate base pairing restraints.
5. projector.dat
6. dummy-contacts-right.dat
7. dummy-contacts-left.dat
8. ss.dat
This file indicates the secondary structure of every residues in the system. In the same order of the input PDB. H: Helix, E: Strand, .: coil. It includes DNA residues all marked as coils.
9. 1a74-BDNA.pdb
This is the input PDB.
