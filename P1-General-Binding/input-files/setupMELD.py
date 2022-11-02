#!/usr/bin/env python
# encoding: utf-8

import numpy as np
from meld.remd import ladder, adaptor, leader
from meld import comm, vault
from meld import system
from meld import parse
from meld.system.restraints import LinearRamp,ConstantRamp
import glob
from restraints import *
import dummy

# number of replicas
N_REPLICAS = 30
# number of steps (units of exchange period)
N_STEPS = 20000
# controles frequence of output
BLOCK_SIZE = 100

with open('projector.dat', 'r') as inputfile:
    projection = inputfile.readlines()
inputfile.close()
projection = [i for i in projection if i != '\n']
num_of_dum = len(projection)


def project(map_list, s):
    '''
    This puts a dummy particle on purine N1 atoms that are involved in contacts.
    '''
    coords = []
    for i in map_list:
        atom = s.index_of_atom(int(i.split()[0]),'N1') - 1
        x,y,z = s.coordinates[atom]
        atom_coords = [x,y,z]
        coords.append(atom_coords)
    
    return coords

def gen_state_templates(index, templates,coor):
    n_templates = len(templates)
    a = system.ProteinMoleculeFromPdbFile(templates[index%n_templates])
    b = system.SystemBuilder(forcefield="ff14sbside")
    dummy_patcher = dummy.DummyPatcher(n_dummies=num_of_dum,coord_dummies=coor)
    c = b.build_system_from_molecules([a],patchers=[dummy_patcher])
    pos = c._coordinates
    vel = np.zeros_like(pos)
    alpha = index / (N_REPLICAS - 1.0)
    energy = 0
    return system.SystemState(pos, vel, alpha, energy,[999,999,999] )

# MAIN routine
def setup_system():
    # create the system starting from coordinates in template.pdb
    templates = glob.glob('*.pdb')
    p = system.ProteinMoleculeFromPdbFile(templates[0])
    b = system.SystemBuilder(forcefield="ff14sbside")
    # load non-standard AA force field params, bonds
    s = b.build_system_from_molecules([p])
    coords = project(projection, s)
    dummy_patcher = dummy.DummyPatcher(n_dummies=num_of_dum,coord_dummies=coords)
    s = b.build_system_from_molecules([p], patchers=[dummy_patcher])
    # Create temperature ladder
    s.temperature_scaler = system.GeometricTemperatureScaler(0.0, 0.3, 300., 450.)
    
    #
    # Preserve protein structure
    #
    const_scaler = s.restraints.create_scaler('constant')
    dist = keep_fixed_distance('protein-contacts.dat',s,scaler=const_scaler)
    s.restraints.add_selectively_active_collection(dist,int(len(dist)))
    #
    # Keep DNA hbonds
    #
    sequenceDNA = readSeq('sequence.dat')
    make_hbond_restraint_file(sequenceDNA,0)
    dist = keep_fixed_distance('hbondsDNA.dat',s,scaler=const_scaler)
    s.restraints.add_selectively_active_collection(dist,int(len(dist)))
    #
    # Keep DNA close to starting conformation
    #
    for i in projection:
        rest = make_cartesian_collections(s, const_scaler, range(int(i.split()[0]),int(i.split()[0])+1),atoms=["N1"], delta=0.6, k=250.)
        s.restraints.add_as_always_active_list(rest)
    #
    # Dummy Contacts
    #
    dummy_cartesian = get_dist_restraints('projector.dat', s, const_scaler, 0, 0.1, 0.1, 0.2, trust_in_group=len(projection))
    s.restraints.add_selectively_active_collection(dummy_cartesian,int(len(dummy_cartesian)))

    dummy_scaler = s.restraints.create_scaler('nonlinear',alpha_min=0.0,alpha_max=1.0, factor=4.0, strength_at_alpha_min=0.5, strength_at_alpha_max=1.0)

    pos1 = s.restraints.create_scaler('linear_positioner',alpha_min=0.0, alpha_max=1.0, pos_min=0.0, pos_max=6.7)
    pos2 = s.restraints.create_scaler('linear_positioner',alpha_min=0.0, alpha_max=1.0, pos_min=0.0, pos_max=6.8)
    pos3 = s.restraints.create_scaler('linear_positioner',alpha_min=0.0, alpha_max=1.0, pos_min=1.5, pos_max=7.2)
    pos4 = s.restraints.create_scaler('linear_positioner',alpha_min=0.0, alpha_max=1.0, pos_min=1.6, pos_max=7.3)

    dummy_contacts_left = get_dist_restraints('dummy-contacts-left.dat', s, dummy_scaler, pos1, pos2, pos3, pos4, trust_in_group=2)
    s.restraints.add_selectively_active_collection(dummy_contacts_left,int(len(dummy_contacts_left)*.5))

    dummy_contacts_right = get_dist_restraints('dummy-contacts-right.dat', s, dummy_scaler, pos1, pos2, pos3, pos4, trust_in_group=2)
    s.restraints.add_selectively_active_collection(dummy_contacts_right,int(len(dummy_contacts_right)*.5))

    #
    # Secondary Structure
    #
    ss_scaler = s.restraints.create_scaler('constant')
    ss_rests = parse.get_secondary_structure_restraints(filename='ss.dat', system=s,ramp=LinearRamp(0,100,0,1), scaler=ss_scaler,
            torsion_force_constant=2.5, distance_force_constant=2.5)
    n_ss_keep = int(len(ss_rests) * 0.96) 
    s.restraints.add_selectively_active_collection(ss_rests, n_ss_keep)

    # Options
    options = system.RunOptions()
    options.implicit_solvent_model = 'gbNeck2'
    options.remove_com = False
    options.use_big_timestep = False # MD timestep (3.3 fs)
    options.use_bigger_timestep = True # MD timestep (4.0 fs)
    options.cutoff = 1.8 # cutoff in nm
    options.soluteDielectric = 1.

    options.use_amap = False # correction to FF12SB
    options.amap_beta_bias = 1.0
    options.timesteps = 11111 # number of MD steps per exchange
    options.minimize_steps = 20000 # init minimization steps

    # create a store
    store = vault.DataStore(s.n_atoms, N_REPLICAS, s.get_pdb_writer(), block_size=BLOCK_SIZE)
    store.initialize(mode='w')
    store.save_system(s)
    store.save_run_options(options)

    # create and store the remd_runner, sets up replica exchange details
    l = ladder.NearestNeighborLadder(n_trials=48)
    policy = adaptor.AdaptationPolicy(2.0, 50, 50)
    a = adaptor.EqualAcceptanceAdaptor(n_replicas=N_REPLICAS, adaptation_policy=policy)
    remd_runner = leader.LeaderReplicaExchangeRunner(N_REPLICAS, max_steps=N_STEPS, ladder=l, adaptor=a)
    store.save_remd_runner(remd_runner)

    # create and store the communicator
    c = comm.MPICommunicator(s.n_atoms, N_REPLICAS)
    store.save_communicator(c)

    # create and save the initial states
    # create and save the initial states, initialize each replica with a different template
    states = [gen_state_templates(i,templates,coords) for i in range(N_REPLICAS)]
    store.save_states(states, 0)

    # save data_store
    store.save_data_store()

    return s.n_atoms

# RUN THE SETUP
setup_system()
