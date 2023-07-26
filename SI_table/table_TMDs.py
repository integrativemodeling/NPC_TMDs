###################################
# Script to summarize all modeling
# information into one table
#
# iecheverria - Salilab - UCSF
# ignacia@salilab.org
###################################

import pandas as pd
import glob
import os
import sys
import numpy as np

sys.path.append('../utils')
from create_summary_table import *
import utils

###########################
# Read files
###########################
modeling_script = '../scripts/mod_TMDs.py'
mmcif_file = '../scripts/IM_TMDs_dimer.cif'

analysis_dir = '../results/'
clustering_dir = os.path.join(analysis_dir,'clustering')
rmf3 = os.path.join(clustering_dir,'cluster.0','cluster_center_model.rmf3')

I = get_input_information(mmcif_file)
input_information = I.get_dictionaries()


R = get_representation(clustering_dir)
representation = R.get_dictionaries()

S = read_modeling_information(modeling_script,
                              analysis_dir,
                              clustering_dir)


sampling = S.get_dictionaries_sampling()


samples = S.get_dictionaries_models()

clustering = S.get_dictionaries_clustering()

#S.update_mmcif_file(mmcif_file)

V = read_validation_information(clustering_dir)
validation = V.get_dictionaries()

#V = read_benchmark_information(clustering_dir)
#benchmark = V.get_dictionaries()

SS = get_software_information(mmcif_file)
software = SS.get_dictionaries()

D = get_data_availability(clustering_dir)
data_availability = D.get_dictionaries()

################################################
# Edit dictionaries
# Entries is dictionaries can be edited to add
# other custom information
################################################
input_information['Experimental data'] = ['14 DSS0']
input_information['Experimental data'].append('Atomic structure prediction from Alphafold2 (AF\_Q12445, AF\_P39685)')
input_information['Prior models'].append('2-fold symmetry derived from cryo-EM structure')

print(representation.keys())

representation['Rigid body (RB) definitions']= ['RB1: Pom152$_{105-130}$', 'RB2: Pom152$_{144-167}$', 'RB3: Pom152$_{176-192}$', 'RB4: Pom152$_{200-212}$', 'RB5: Pom34$_{44-86}$', 'RB6: Pom34$_{89-110}$']

representation['Resolution of structured components'] = '1 [R1] residue per bead'
representation['Resolution of disordered regions'] = '2 [R2] residues per bead'

representation['Composition (number of copies)'] = 4 


representation['Spatial restraints encoded into scoring function'] = representation.pop('Spatial restraints encoded into scoring function')

representation['Spatial restraints encoded into scoring function'].append('Cross-link restraints; applied to the R1 representation')
#representation['Spatial restraints encoded into scoring function'].append('Sequence connectivity and excluded volume')



sampling['CPU time'] = ['22 hours on 80 processors']


validation['Percent of sequence connectivity restraints satisfied per structure'] = ['99 \%']
validation['Percent cross-link restraints satisfied by ensemble'] = ['99 \%']
validation['Percent of excluded volume restraints satisfied per structure'] = ['99 \%']



print(samples)

#samples['Number of models after equilibration']= '300000/300000'
#samples['Number of models that satisfy the input information']= '92281/51183'
#samples['Number of structures in samples A/B']= '52831,39450/18722,32461' 
#samples['p-value of non-parametric Kolmogorov-Smirnov two-sample test']= '0.446/0.34'
#samples['Kolmogorov-Smirnov two-sample test statistic, D']= '0.0,1E-16'

#samples['Number of models after equilibration (isolated/in situ)']=samples.pop('Number of models after equilibration')
#samples['Number of models that satisfy the input information (isolated/in situ)']=samples.pop('Number of models that satisfy the input information')
#samples['Number of structures in samples A,B (isolated/in situ)']=samples.pop('Number of structures in samples A/B')
#samples['p-value of non-parametric Kolmogorov-Smirnov two-sample test (isolated/in situ)']=samples.pop('p-value of non-parametric Kolmogorov-Smirnov two-sample test')
#samples['Kolmogorov-Smirnov two-sample test statistic (D) (isolated/in situ)']=samples.pop('Kolmogorov-Smirnov two-sample test statistic, D')


###################
print('-------------')
print(clustering)
#clustering['Cluster precisions']=['cluster 1 : 87.0 $\\%$']
#clustering['Cluster populations']=clustering.pop('Cluster populations')

sampling['Replica exchange temperature range'] =  ['1.0 - 4.0']

#clustering['Sampling precision']='14.66/12.6 \\AA'
#clustering["Homogeneity of proportions $\\chi^2$ test (p-value)/Cramer’s V value"]="1.000 (0.000)/ 1.000 (0.000) (thresholds: p-value$>$0.05 OR Cramer's V$<$0.1)"
#"1.000 (0.000)/ 1.000 (0.000)(thresholds: p-value$>$0.05 OR Cramer's V$<$0.1)"]
#clustering['Number of clusters']='1/1'
#clustering['Cluster populations']='Cluster 1: 99/98 \%'
#clustering['Cluster precisions']='Cluster 1: 15.42/10.1 \\AA'
#clustering['Average cross-correlation between localization probability densities of samples A and B']='Cluster 1: 0.84/0.92'

#clustering['Sampling precision']=clustering.pop('Sampling precision')
#clustering['Homogeneity of proportions $\\chi^2$ test p-value (Cramer’s V value)']=clustering.pop('Homogeneity of proportions $\\chi^2$ test (p-value)/Cramer’s V value')
#clustering['Number of clusters']=clustering.pop('Number of clusters')
#clustering['Cluster populations']=clustering.pop('Cluster populations')
#clustering['Cluster precisions']=clustering.pop('Cluster precisions')
#clustering['Average cross-correlation between localization probability densities of samples A and B']=clustering.pop('Average cross-correlation between localization probability densities of samples A and B')


software['Modeling scripts'] = ['https://github.com/integrativemodeling/NPC\_TMD']
software['Structure prediction'] = ['AlphaFold2']
software['Visualization and plotting'] = ['UCSF Chimera']

################################################
# Convert ordered dictionaries 
# into lists
################################################
input_information_list = dict_to_list(input_information)
representation_list = dict_to_list(representation)
sampling_list = dict_to_list(sampling)
samples_list = dict_to_list(samples)
clustering_list = dict_to_list(clustering)
validation_list = dict_to_list(validation)
software_list = dict_to_list(software)
data_availability_list = dict_to_list(data_availability)

print(sampling_list)


################################################
# Compile all information
# 
################################################
variable_dict = {'complex': 'Pom34-Pom152 TMDs',
                 'number':10,
                 'input_information': input_information_list, 
                 'representation': representation_list,
                 'sampling': sampling_list,
                 'samples': samples_list,
                 'clustering':clustering_list,
                 'validation':validation_list,
                 'software':software_list,
                 'data':data_availability_list}

################################################
# Generate tex, pdf file
################################################
template = utils.get_template('../utils/SI_template.tex')
utils.compile_pdf_from_template(template, variable_dict, './table_SI_TMDs.pdf')

exit()
