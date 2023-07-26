# Integrative structure determination of the Pks13 dimer

## Summary

A structural model of the Pks13 dimer was computed by integrative modeling  based on the 57 DSSO chemical cross-links and structural models of the components of the Pks13 dimer. These components include the atomic structural model based on the cryo-EM maps (EMD-26574), de novo AlphaFold2 predictions of domains unresolved in the cryo-EM structure and flexible linker regions. A model of the Pks13 dimer was computed by satisfying this input information to the best possible degree using IMP. The resulting ensemble of acceptable models satisfies 91% of the cross-links. 

## List of files and directories:

- `data` All data used for integrative modeling, including the cross-links, the partial Pks13 structure obtained using cryo-EM, and the AlphaFold model. 
- `scripts` PMI modeling script (`mod_symmetry.py`) to model the symmetry full-lenght Pks13 dimer.
- `analysis` Scripts to analyze the simulations 
- `results` All the relevant results from integrative modeling, including the distance statistics for the cross-links and clustering of the results.
- `SI_table` Scripts to generate a table summarizing the integrative modeling protocols.
- `utils` Template and code to generate the <em>Supporting information</em> table summarizing the integrative modeling protocol.

## Information

*Author (s)*: Ignacia Echeverria

_License_: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/) This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License.

_Publications_: Structure and dynamics of the essential endogenous mycobacterial polyketide synthase Pks13. Sun Kyung Kim, Miles Sasha Dickinson,, Janet Finer-Moore, Ziqiang Guan, Robyn M. Kaake, Ignacia Echeverria, Jen Chen, Ernst H. Pulido, Andrej Sali, Nevan J. Krogan, Oren S. Rosenberg, Robert M. Stroud, Nature Structural & Molecular Biology. 2022