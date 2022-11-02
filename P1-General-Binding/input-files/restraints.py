import numpy as np
from meld.remd import ladder, adaptor, leader
from meld import comm, vault
from meld import system
from meld import parse
from meld.system.restraints import LinearRamp,ConstantRamp


def make_cartesian_collections(s, scaler, residues, atoms = ['CA'], delta=0.35, k=250.):
    '''
    Used to restrain atoms to their initial coordinates.
    '''
    cart = []
    atoms_restraint = atoms
    sequence = [(i,j) for i,j in zip(s.residue_numbers,s.residue_names)]
    sequence = sorted(set(sequence))
    sequence = dict(sequence)
    for i in residues:
        if sequence[i] in 'ACE':
            atoms_restraint = ['CH3','C','O']
        elif sequence[i] in 'NME':
            atoms_restraint = ['CH3','N']
        else:
            atoms_restraint = atoms
        for b in atoms_restraint:
            try:
                atom_index = s.index_of_atom(i,b) - 1
                x,y,z = s.coordinates[atom_index]/10.
                rest = s.restraints.create_restraint('cartesian',scaler, res_index=i, atom_name=b,
                    x=x, y=y, z=z, delta=delta,force_const=k)
                cart.append(rest)
            except:
                pass
    return cart


def keep_fixed_distance(filename, s, scaler):
    '''
    Used to keep a specific contact between residue 1 atom1 residue2 atom2 at a distance that is +- 1 the original distance
    '''
    dists = []
    rest_group = []
    lines = open(filename).read().splitlines()
    lines = [line.strip() for line in lines]
    for line in lines:
        if not line:
            dists.append(s.restraints.create_restraint_group(rest_group, 1))
            rest_group = []
        else:
            cols = line.split()
            i = int(cols[0])
            name_i = cols[1]
            j = int(cols[2])
            name_j = cols[3]
            dist = float(cols[4]) / 10.

            rest = s.restraints.create_restraint('distance', scaler,LinearRamp(0,100,0,1),
                                                 r1=dist-0.2, r2=dist-0.1, r3=dist+0.1, r4=dist+0.2, k=250,
                                                 atom_1_res_index=i, atom_2_res_index=j,
                                                 atom_1_name=name_i, atom_2_name=name_j)
            rest_group.append(rest)
    return dists

def make_hbond_restraint_file(sequence,offset,name='hbondsDNA.dat'):
    '''
    Used to restrain DNA base pairing.
    '''
    offset = int(offset)
    tot = len(sequence)
    print(sequence,tot)

    AT = '''{0} N1 {1} N3 3.0
    
{0} N6 {1} O4 3.0
    
'''
    
    CG = '''{0} O2 {1} N2 3.0
    
{0} N3 {1} N1 3.0
    
{0} N4 {1} O6 3.0
    
'''
    
    GC = '''{0} N2 {1} O2 3.0
    
{0} N1 {1} N3 3.0
    
{0} O6 {1} N4 3.0
    
'''
    
    TA = '''{0} N3 {1} N1 3.0
    
{0} O4 {1} N6 3.0
    
'''
    
    values = {'A':AT,'C':CG,'G':GC,'T':TA}
    
    with open(name,'w') as fo:
        for i,j in enumerate(sequence):
            ni = i + 1
            ncomplementary = tot*2 - ni +1
            ni = ni + offset
            ncomplementary += offset
            fo.write(values[j].format(ni,ncomplementary))


def readSeq(fname):
    '''
    Reads a Levitt like sequence file with the following format
    SEQ A >CCTTGGCTGACGTCAGCCAAG
    SEQ B >CTTGGCTGACGTCAGCCAAGG
    We want the first strand

    Or read a sequence file CCCCTACACT
    '''
    fi = open(fname,'r').readlines()[0].rstrip()
    if 'SEQ' in fi:
        return fi.split(">")[1]
    else:
        return fi


def get_dist_restraints(filename, s, scaler, p1, p2, p3, p4, trust_in_group=2):
    dists = []
    rest_group = []
    lines = open(filename).read().splitlines()
    lines = [line.strip() for line in lines]
    for line in lines:
        if not line:
            dists.append(s.restraints.create_restraint_group(rest_group, trust_in_group))
            rest_group = []
        else:
            cols = line.split()
            i = int(cols[0])
            name_i = cols[1]
            j = int(cols[2])
            name_j = cols[3]
            # dist = float(cols[4]) / 10.

            rest = s.restraints.create_restraint('distance', scaler,LinearRamp(0,100,0,1),
                                                 r1=p1, r2=p2, r3=p3, r4=p4, k=250,
                                                 atom_1_res_index=i, atom_2_res_index=j,
                                                 atom_1_name=name_i, atom_2_name=name_j)
            rest_group.append(rest)
    return dists
