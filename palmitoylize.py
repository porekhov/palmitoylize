import sys, os, optparse, re

parser = optparse.OptionParser()
parser.add_option('-p','--pdb', dest='pdb', help = 'PDB file')
parser.add_option('-i','--itp', dest='itp', help = 'itp file')
parser.add_option('-d','--dir', dest='direct', help = 'direction of palmitoyl tail with respect to the membrane normal: 1 or -1')
parser.add_option('-c','--cys', dest='cys', help = 'list of cys res indices to modify')

(options,args) = parser.parse_args()
pdb_file = os.path.abspath(options.pdb)
itp_file = os.path.abspath(options.itp)
direction = int(options.direct)
cys2cyp = str(options.cys)

class Atom:
    '''Atom object'''
    def __init__(self, id, type, resid, resname, name, cgnr, charge):
        '''Init function for Atom object'''
        self.id = int(id)
        self.resid = int(resid)
        self.resname = str(resname)
        self.type = str(type)
        self.name = str(name)
        self.charge = float(charge)
        self.cgnr = int(cgnr)
    
    def fmt(self):
        '''formatted line for Atom object'''
        return "%3d %5s %3d %5s %5s %3d %1.3f" \
            % (self.id, self.type, self.resid, self.resname, self.name, self.cgnr, self.charge)

class AtomPDB:
    """This class represents atom properties"""
    def __init__(self, xyz = [], atype = '', chain_id = 'A', res_type = '', res_num = 0, num = 0):
        '''Init function for Atom object'''
        self.xyz = xyz
        self.type = atype
        self.chain_id = chain_id
        self.res_type = res_type
        self.res_num = res_num
        self.num = num
        
    def pdb_str(self):
        '''Print PDB formatted line for Atom object'''
        return "%6s%5s %4s %3s %1s%4d    %8.3f%8.3f%8.3f" \
            % ("ATOM  ", self.num, self.type, self.res_type, self.chain_id, self.res_num, self.xyz[0], self.xyz[1], self.xyz[2])
               
def AtomFromPdbLine(line):
    """returns an Atom object from an atom line in a pdb file."""
    atom = AtomPDB()
    atom.num = int(line[6:11])
    atom.type = line[12:16].strip(" ")
    atom.res_type = line[17:20]
    atom.chain_id = line[21]
    atom.res_num = int(line[22:26])
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    atom.xyz = [x, y, z]
    return atom

class MoleculePDB:
    '''This class serves as a container for atom objects and reads/writes pdb files'''
    def __init__(self, pdb=""):
        '''Init function for Molecule object'''
        self.id = ''
        self._atoms = []
        if pdb:
            self.read_pdb(pdb)

    def atom(self, i):
        '''Return atom i'''
        return self._atoms[i]
    
    def clear(self):
        '''Clear all Atoms from Molecule'''
        for atom in self._atoms:
            del atom
        self._atoms = []

    def insert_atom(self, atom):
        '''Add Atom to Molecule'''
        self._atoms.append(atom)
  
    def read_pdb(self, fname):
        '''Load Atoms from pdb file'''
        self.clear()
        for line in open(fname, 'r').readlines():
            if line.startswith("ATOM"):
                atom = AtomFromPdbLine(line);
                if len(self._atoms) == 1:
                    self.id = atom.chain_id
                self.insert_atom(atom)

atoms, cysteines = [], []
cys2cyp = [i.strip() for i in cys2cyp.split(',')]
tot_atoms = 0
cur_sec = '' # currently processed section of an itp file
        
for line in open(itp_file, 'r').readlines():
    if line.replace(" ", "")[0] != ';' and line.replace(" ", "")[0] != "#" \
                                               and re.search('[a-zA-Z0-9]', line):
        fields = line.split()
        if line.replace(" ", "").startswith('[moleculetype]'):
            cur_section = 'mt'
            continue
        elif line.replace(" ", "").startswith('[atoms]'):
            cur_section = 'at'
            continue
        elif line.replace(" ", "").startswith('[bonds]'):
            cur_section = 'b'
            continue
        elif line.replace(" ", "").startswith('[angles]'):
            cur_section = 'an'
            continue
        elif line.replace(" ", "").startswith('[dihedrals]'):
            cur_section = 'd'
            continue
        elif line.replace(" ", "").startswith('[constraints]'):
            cur_section = 'c'
            continue
        elif line.replace(" ", "").startswith('[exclusions]'):
            cur_section = 'e'
            continue
        elif line.replace(" ", "").startswith('[position_restraints]'):
            cur_section = 'pc'
            continue
        elif line.replace(" ", "").startswith('[virtual_sites2]'):
            cur_section = 'vs2'
            continue
            # handling atoms
        if cur_section == 'at' and fields[3] == 'CYS' and fields[4] == 'SC1' and fields[2] in cys2cyp:
            cysteines.append(Atom(fields[0], fields[1], fields[2], 'CYP', \
                                                fields[4], fields[5], fields[6]))
            atoms.append(Atom(fields[0], fields[1], fields[2], 'CYP', \
                                                fields[4], fields[5], fields[6]))
            tot_atoms += 1
        elif cur_section == 'at' and fields[3] == 'CYS' and fields[4] == 'BB' and fields[2] in cys2cyp:
            atoms.append(Atom(fields[0], fields[1], fields[2], 'CYP', \
                                                fields[4], fields[5], fields[6]))
            tot_atoms += 1
        elif cur_section == 'at':
            atoms.append(Atom(fields[0], fields[1], fields[2], fields[3], \
                                                fields[4], fields[5], fields[6]))
            tot_atoms += 1

extr_at_count = 1
atom_names = ['C5', 'C4', 'C3', 'C2']
atom_displace = [3.9, 8.6, 13.3, 17.4]
new_atoms, new_pdb, new_bonds, new_angles = [], [], [], []
pdb = MoleculePDB(pdb_file)

for cys in cysteines:
    added_atoms = []
    for i in range(4):
        atom_id = tot_atoms + extr_at_count
        added_atoms.append(atom_id)
        new_atoms.append("%3d %5s %3d %5s %5s %3d %1.3f" % (atom_id, 'C1', cys.resid, cys.resname, atom_names[i], atom_id, 0))
        cys_x = pdb.atom(cys.id - 1).xyz[0]
        cys_y = pdb.atom(cys.id - 1).xyz[1]
        cys_z = pdb.atom(cys.id - 1).xyz[2]
        atom_to_add = AtomPDB([cys_x, cys_y, cys_z + atom_displace[i]*direction], atom_names[i], pdb.atom(cys.id - 1).chain_id, 'CYP', cys.resid, atom_id)
        new_pdb.append(atom_to_add.pdb_str())
        extr_at_count += 1
        

    new_bonds.append("%3d %3d %1d %5.3f %5.3f" % (cys.id, added_atoms[0], 1, 0.390, 5000))
    new_bonds.append("%3d %3d %1d %5.3f %5.3f" % (added_atoms[0], added_atoms[1], 1, 0.470, 4500))
    new_bonds.append("%3d %3d %1d %5.3f %5.3f" % (added_atoms[1], added_atoms[2], 1, 0.470, 4500))
    new_bonds.append("%3d %3d %1d %5.3f %5.3f" % (added_atoms[2], added_atoms[3], 1, 0.410, 5500))

    new_angles.append("%3d %3d %3d %1d %5.3f %5.3f" % (cys.id - 1, cys.id, added_atoms[0], 2, 150.000, 35.0))
    new_angles.append("%3d %3d %3d %1d %5.3f %5.3f" % (cys.id, added_atoms[0], added_atoms[1], 2, 125.000, 25.0))
    new_angles.append("%3d %3d %3d %1d %5.3f %5.3f" % (added_atoms[0], added_atoms[1], added_atoms[2], 2, 180.000, 25.0))
    new_angles.append("%3d %3d %3d %1d %5.3f %5.3f" % (added_atoms[1], added_atoms[2], added_atoms[3], 2, 180.000, 28.0))
    
print('[ atoms ]')
for atom in atoms:
    print(atom.fmt())
print('; cys-palmitoyl beads')
print('\n'.join(new_atoms))
print('\nAdd to the [ bonds ] section:\n')
print('\n'.join(new_bonds))
print('\nAdd to the [ angles ] section:\n')
print('\n'.join(new_angles))
print('\nAdd to PDB file\n')
print('\n'.join(new_pdb))
