# palmitoylize
a script for palmitoylation of coarse-grained MARTINI protein models

**Usage:**
palmitoylize -p protein.pdb -i protein.itp -d '1 or -1' -c 'cysteine residue indices to modify as they appear in protein.itp comma separated'

protein.itp can be obtained using martinize.py script

protein.pdb should contain a membrane protein aligned such that Z direction corresponds to the membrane normal. If palmytoil tails should expand in the positive direction of Z use -d 1 otherwise use -d -1.
