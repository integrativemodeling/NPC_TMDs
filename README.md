# Integrative structure determination of the Pom34-152 transmembrane

## Summary

A structural model of the Pom34-152 transmembrane attachment site dimer was computed by integrative modeling  based on the cryo-EM density maps (EMD-41117) and the AlphaFold2 structure predcition of the monomers. A model of the dimer Pom34-152 transmembrane domains (TMDs) was computed by satisfying this input information to the best possible degree using IMP. 

## List of files and directories:

- `data` All data used for integrative modeling, including the cross-links, EM density map, and the AlphaFold2 models. 
- `scripts` PMI modeling script (`mod_TMDs.py`) to model the symmetric Pom34-152 dimer.
- `analysis` Scripts to analyze the simulations 
- `results` All the relevant results from integrative modeling, including the distance statistics for the cross-links and clustering of the results.
- `SI_table` Scripts to generate a table summarizing the integrative modeling protocols.
- `utils` Template and code to generate the <em>Supporting information</em> table summarizing the integrative modeling protocol.

## Information

*Author (s)*: Ignacia Echeverria

_License_: [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/) This work is licensed under the Creative Commons Attribution-ShareAlike 4.0 International License.

_Publications_: Implications of a multiscale structure of the yeast Nuclear Pore Complex. Akey CA, Echeverria I, Ouch C, Nudelman I, Shi Y, Wang J, Weiss TM, Shi Y, Chait BT, Sali A,Fernandez-Martinez J, Rout MP. Molecular Cell. 2023
