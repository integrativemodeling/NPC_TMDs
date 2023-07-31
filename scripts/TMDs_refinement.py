import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
import IMP.pmi.restraints.em 
import IMP.pmi.restraints.basic
import IMP.pmi.io.crosslink
import IMP.pmi.restraints.npc
import IMP.bayesianem
import IMP.bayesianem.restraint
import math
import time
import numpy as np
import os

from npc import *

top_dir = '/Users/iecheverria/Dropbox/UCSF/yeast_npc/modeling/mod_franken/mod_pomring/'
top_dir = '/wynton/home/sali/ignacia/NPC/modeling_2023/TMDs/'

# Create System and State
mdl = IMP.Model()
s = IMP.pmi.topology.System(mdl)
st = s.create_state()

mols = []
domains = {'Pom152':[[105,130], [144,167], [176, 192], [200, 212]],
           'Pom34':[[44,86],[89,110],[122,150],[222,237]]}

pdb = {'Pom152':'data/selected.rebuilt.pdb',
       'Pom34':'data/selected.rebuilt.pdb'}

chain_id = ['A','D']
n_densities = 10
colors = ['gray','red' ]

for k, (prot, doms) in enumerate(domains.items()):
    seqs = IMP.pmi.topology.Sequences(top_dir+f'data/{prot}.fasta')
    mol = st.create_molecule(prot,sequence=seqs[prot],chain_id=chain_id[k])
    pdb_file = os.path.join(top_dir,pdb[prot])
    struct = []
    for i, d in enumerate(doms):
        print('-----', d)
        a1 = mol.add_structure(pdb_file,
                           chain_id=chain_id[k],
                           res_range = d)

        # Add structured part representation and then build
        mol.add_representation(a1,
                               density_residues_per_component=n_densities,
                               density_voxel_size=3.0,
                               resolutions=[1],
                               density_prefix = os.path.join(f'{top_dir}/data/em_data',f'{prot}_{i}'),
                               color = colors[k])
        struct += a1
    mol.add_representation(mol[0:250]-struct,resolutions=[2], color = colors[k])
    mols.append(mol)
        
# Clone
clones = []
chains='DE'
for i,mol in enumerate(mols):
    clone = mol.create_clone(chains[i])
    clones.append(clone)

mols = mols + clones
    
hier = s.build()

mdl.update() # propagates coordinates

out = IMP.pmi.output.Output()
out.init_rmf("sym_ini.rmf3", [hier])
out.write_rmf("sym_ini.rmf3")
out.close_rmf("sym_ini.rmf3")
############################################
# DOFs
############################################             
dof = IMP.pmi.dof.DegreesOfFreedom(mdl)

for i, mol in enumerate(mols):
    dof.create_flexible_beads(mol.get_non_atomic_residues())
    copy_index =  IMP.atom.Copy(mol.get_hierarchy()).get_copy_index()
    print(mol.get_name(), copy_index)
    for dom in domains[mol.get_name()]:
        sel = IMP.atom.Selection(hier,
                                 molecule=mol.get_name(),
                                 residue_indexes=range(dom[0],dom[1]+1,1),
                                 copy_index=copy_index,
                                 resolution=IMP.atom.ALL_RESOLUTIONS).get_selected_particles()
        sel_densities = IMP.atom.Selection(hier,
                                           molecule=mol.get_name(),
                                           residue_indexes=range(dom[0],dom[1]+1,1),
                                           copy_index=copy_index,
                                           representation_type=IMP.atom.DENSITIES).get_selected_particles()

        dof.create_rigid_body(sel+sel_densities)
   
print(dof.get_movers())

# Symmetry

center = IMP.algebra.Vector3D([0,0,0])

# Axial symmetry
rot = IMP.algebra.get_rotation_about_axis([0,1,0],math.pi)
transform = IMP.algebra.get_rotation_about_point(center,rot)
dof.constrain_symmetry(mols[0],mols[2],transform)
dof.constrain_symmetry(mols[1],mols[3],transform)

#############################
# Connectivity restraint
#############################
output_objects = []

for n, mol in enumerate(mols):
    cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
    cr.set_label(mol.get_name()+'.'+str(n))
    cr.add_to_model()
    output_objects.append(cr)

##########################
# Membrane binding
##########################
#tor_th      = 45.0
#tor_th_ALPS = 12.0
#tor_R       = 390.0 + 225.0
#tor_r       = 150.0 - tor_th/2.0
#tor_r_ALPS  = 150.0 - tor_th_ALPS/2.0
#msl_sigma   = 1.0
#msl_weight  = 10.0

tor_th      = 45.0
tor_th_ALPS = 12.0
tor_R       = 390.0 + 180.0
tor_r       = 150.0 - tor_th/2.0
tor_r_ALPS  = 150.0 - tor_th_ALPS/2.0
msl_sigma   = 2.0
msl_weight  = 1.0
###########################
# Membrane binding
###########################

transmembrane_sel = [(106,133,'Pom152',0),
                     (144,167,'Pom152',0),
                     (176,194,'Pom152',0),
                     (57,87,'Pom34',0),
                     (122,150,'Pom34',0)]

for sel in transmembrane_sel:
    msl = IMP.pmi.restraints.npc.MembraneSurfaceLocationRestraint(hier=hier,
                                                                  protein=sel,
                                                                  tor_R=tor_R,
                                                                  tor_r=tor_r,
                                                                  tor_th=tor_th,
                                                                  sigma=msl_sigma,
                                                                  resolution = 1)
    msl.set_weight(msl_weight)
    msl.add_to_model()
    #msl.set_label('%s.%s'%(sel[2],sel[0]))
    output_objects.append(msl)
    print('Membrane binding restraint:', msl.evaluate())

perinuclear_sel = [(134,143,'Pom152',0),
                   (195,210,'Pom152',0),
                   (88,89,'Pom34',0),
                   (111,118,'Pom34',0)]
                      
    
for sel in perinuclear_sel:
    print('Applying membrane localization restraint:', sel)
    msl_1 = PerinuclearVolumeLocationRestraint(hier,
                                                                      protein=sel,
                                                                      tor_R=tor_R,
                                                                      tor_r=tor_r,
                                                                      tor_th=tor_th,
                                                                      sigma=msl_sigma,
                                                                      resolution = 10,
                                                                      label=f'{sel[2]}.{sel[0]}')
    msl_1.set_weight(msl_weight)
    msl_1.add_to_model()
    output_objects.append(msl_1)
    print('Membrane binding restraint ready', msl_1.get_output())

poreside_sel = [(1,43,'Pom34',0),
                (158,299,'Pom34',0),
                (1,104,'Pom152',0),
                (167,175,'Pom152',0)]
    
    
for sel in poreside_sel:
    print('Applying membrane localization restraint:', sel)
    msl_1 = PoreSideVolumeLocationRestraint(hier,
                                                                   protein=sel,
                                                                   tor_R=tor_R,
                                                                   tor_r=tor_r,
                                                                   tor_th=tor_th,
                                                                   sigma=msl_sigma,
                                                                   resolution = 10,
                                                                   label=f'{sel[2]}.{sel[0]}')
    msl_1.set_weight(msl_weight)
    msl_1.add_to_model()
    output_objects.append(msl_1)
    print('Membrane binding restraint ready', msl_1.get_output())

############################################
# EM
############################################     

densities = IMP.atom.Selection(hier,
                               representation_type=IMP.atom.DENSITIES).get_selected_particles()
    
IMP.isd.gmm_tools.write_gmm_to_map(densities, "test.mrc", voxel_size=3.0, fast=True)

print('densities', densities, len(densities))
gem = IMP.bayesianem.restraint.GaussianEMRestraintWrapper(densities,
                                                          top_dir+'data/em_data/run27h-c2-50-300bf_locres_filt-mask2-2-norm-100-zone16-TMDs-dust_dsfact1_ng60_oriented.txt',
                                                          scale_target_to_mass=True)
gem.set_label("EM_membrane")
#gem.add_target_density_to_hierarchy(st)
gem.add_to_model()
gem.set_weight(20.0)
output_objects.append(gem)
#gem.center_model_on_target_density(st)

t0 = gem.evaluate()

############################
# Excluded Volume
############################

# Create excluded volume for all particles
evr1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=mols)
evr1.add_to_model()
evr1.set_weight(1.0)
#evr1.set_label('intra')
output_objects.append(evr1)

##########################
# Inter-molecular XLs
##########################
resi_XLs = [62, 301, 351]
resi_XLs = [(62,'Pom152'), (44,'Pom34')]
w_intra_xls = 5.
include_intra_XLs = True
if include_intra_XLs:
    for (r, prot) in resi_XLs:
        dist_min = 10.0
        dist_max = 30.0
        ixl = IMP.pmi.restraints.basic.DistanceRestraint(root_hier = hier,
                                                         tuple_selection1=(r,r,prot,0),
                                                         tuple_selection2=(r,r,prot,1),
                                                         distancemin=dist_min,
                                                         distancemax=dist_max,
                                                         label=f'XLs_inter_{r}')
        ixl.set_weight(w_intra_xls)
        ixl.add_to_model()
        output_objects.append(ixl)
        print('Intra molecular XLs:', ixl.get_output())

###########################
# Chemical crosslinks
###########################
# INITIALIZE DB
rmf_restraints = []
include_XLs = True
w_xls = 5.
if include_XLs:
    cldbkc=IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
    cldbkc.set_protein1_key("Protein1")
    cldbkc.set_protein2_key("Protein2")
    cldbkc.set_residue1_key("Residue1")
    cldbkc.set_residue2_key("Residue2")
    #cldbkc.set_unique_id_key("UniqueID")
    #cldbkc.set_psi_key("Psi")

    # XLs RESTRAINT
    cldb=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
    cldb.create_set_from_file(top_dir+"data/XLs_all_2020.csv")

    xl1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=hier,
                                                                                database=cldb,
                                                                                resolution=1,
                                                                                length=21.0,
                                                                                slope=0.01)
    xl1.add_to_model()
    xl1.set_weight(w_xls)
    rmf_restraints.append(xl1)
    output_objects.append(xl1)
        
###########################
# Randomize configurations
###########################

sel = IMP.atom.Selection(hier).get_selected_particles()

#IMP.pmi.tools.shuffle_configuration(sel,
#                                    bounding_box=((-150, -540, -20), (150, -440, 20)),
#                                    avoidcollision_rb=False)

mdl.update()
IMP.isd.gmm_tools.write_gmm_to_map(densities, "test_shuffle.mrc", voxel_size=3.0, fast=True)

############################
# Sampling
############################

mc1 = IMP.pmi.macros.ReplicaExchange(mdl,
                                      root_hier=hier,                       
                                      #crosslink_restraints=rmf_restraints,       
                                      monte_carlo_sample_objects=dof.get_movers(),  
                                      global_output_directory="output/",
                                      output_objects=output_objects,
                                      replica_exchange_maximum_temperature=4.0,
                                      monte_carlo_steps=20,
                                      number_of_frames=30000,
                                      number_of_best_scoring_models=0)

mc1.execute_macro()
rex1 = mc1.get_replica_exchange_object()

exit()
