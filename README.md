# palmitoylize
a script for palmitoylation of coarse-grained MARTINI protein models

**Usage:**
palmitoylize -p protein.pdb -i protein.itp -d '1 or -1' -c 'cysteine residue indices to modify as they appear in protein.itp comma separated'

protein.itp can be obtained using martinize.py script.

protein.pdb should contain a coarse-grained membrane protein (no water/ions/lipids) aligned such that Z direction corresponds to the membrane normal. If palmytoil tails should expand in the positive direction of Z use -d 1 otherwise use -d -1. This file can be also obtained by martinize.py.

The script returns the complete [ atoms ] section of topolgy with CYS names/bead types modified accoringly and the extra atoms added in the end of this section. This section along with the returned additional bond and agnle terms should be manually added to the original itp file.

Finally, the additional beads printed in the pdb format should be added at the end of the original pdb file.
