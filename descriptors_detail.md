# Descriptors Used in Water Solubility Prediction Study

This document provides a detailed list of all descriptors used in the water solubility prediction study. Each descriptor is briefly described to give context on its relevance and role in the model.

## Table of Contents
1. [Basic Descriptors](#basic-descriptors)
2. [Functional Group Descriptors](#functional-group-descriptors)
3. [Feature-Engineered Descriptors](#feature-engineered-descriptors)
4. [Descriptors Used by Sorkun](#descriptors-used-by-sorkun)

---

## Basic Descriptors
1. **MaxAbsEStateIndex**: The maximum absolute value of the E-state index, a measure of electron distribution.
2. **MaxEStateIndex**: The maximum E-state index in the molecule.
3. **MinAbsEStateIndex**: The minimum absolute value of the E-state index.
4. **MinEStateIndex**: The minimum E-state index in the molecule.
5. **qed**: Quantitative estimate of drug-likeness, a metric used to assess the drug-likeness of a molecule.
6. **SPS**: Synthetic Accessibility Score, predicting the ease of synthesis for a compound.
7. **MolWt**: Molecular Weight, the total mass of the molecule.
8. **HeavyAtomMolWt**: The molecular weight of heavy atoms (non-hydrogen atoms).
9. **ExactMolWt**: The exact molecular weight considering isotopic distribution.
10. **NumValenceElectrons**: Total number of valence electrons in the molecule.
11. **NumRadicalElectrons**: Number of radical electrons in the molecule.
12. **MaxPartialCharge**: Maximum partial charge on any atom in the molecule.
13. **MinPartialCharge**: Minimum partial charge on any atom in the molecule.
14. **MaxAbsPartialCharge**: Maximum absolute value of the partial charge on any atom.
15. **MinAbsPartialCharge**: Minimum absolute value of the partial charge on any atom.
16. **FpDensityMorgan1**: Morgan fingerprint density with radius 1, a circular fingerprint used to encode molecular structure.
17. **FpDensityMorgan2**: Morgan fingerprint density with radius 2.
18. **FpDensityMorgan3**: Morgan fingerprint density with radius 3.
19. **BCUT2D_MWHI**: BCUT metric using molecular weight as the property.
20. **BCUT2D_MWLOW**: BCUT metric for low molecular weight property.
21. **BCUT2D_CHGHI**: BCUT metric using high charge as the property.
22. **BCUT2D_CHGLO**: BCUT metric for low charge property.
23. **BCUT2D_LOGPHI**: BCUT metric using high LogP as the property.
24. **BCUT2D_LOGPLOW**: BCUT metric for low LogP property.
25. **BCUT2D_MRHI**: BCUT metric using high molar refractivity as the property.
26. **BCUT2D_MRLOW**: BCUT metric for low molar refractivity property.
27. **AvgIpc**: Average Information Content (Ipc) of the molecule.
28. **BalabanJ**: Balaban’s J index, a topological descriptor.
29. **BertzCT**: Bertz complexity index, a measure of molecular complexity.
30. **Chi0**: The first-order connectivity index.
31. **Chi0n**: The first-order connectivity index with nitrogen.
32. **Chi0v**: The first-order valence connectivity index.
33. **Chi1**: The second-order connectivity index.
34. **Chi1n**: The second-order connectivity index with nitrogen.
35. **Chi1v**: The second-order valence connectivity index.
36. **Chi2n**: The third-order connectivity index with nitrogen.
37. **Chi2v**: The third-order valence connectivity index.
38. **Chi3n**: The fourth-order connectivity index with nitrogen.
39. **Chi3v**: The fourth-order valence connectivity index.
40. **Chi4n**: The fifth-order connectivity index with nitrogen.
41. **Chi4v**: The fifth-order valence connectivity index.
42. **HallKierAlpha**: Hall-Kier alpha modification of the shape index.
43. **Ipc**: Information content index, a molecular descriptor.
44. **Kappa1**: First Kappa shape index.
45. **Kappa2**: Second Kappa shape index.
46. **Kappa3**: Third Kappa shape index.
47. **LabuteASA**: Labute’s approximation to molecular surface area.
48. **PEOE_VSA1**: PEOE charge on molecular surface area, bin 1.
49. **PEOE_VSA10**: PEOE charge on molecular surface area, bin 10.
50. **PEOE_VSA11**: PEOE charge on molecular surface area, bin 11.
51. **PEOE_VSA12**: PEOE charge on molecular surface area, bin 12.
52. **PEOE_VSA13**: PEOE charge on molecular surface area, bin 13.
53. **PEOE_VSA14**: PEOE charge on molecular surface area, bin 14.
54. **PEOE_VSA2**: PEOE charge on molecular surface area, bin 2.
55. **PEOE_VSA3**: PEOE charge on molecular surface area, bin 3.
56. **PEOE_VSA4**: PEOE charge on molecular surface area, bin 4.
57. **PEOE_VSA5**: PEOE charge on molecular surface area, bin 5.
58. **PEOE_VSA6**: PEOE charge on molecular surface area, bin 6.
59. **PEOE_VSA7**: PEOE charge on molecular surface area, bin 7.
60. **PEOE_VSA8**: PEOE charge on molecular surface area, bin 8.
61. **PEOE_VSA9**: PEOE charge on molecular surface area, bin 9.
62. **SMR_VSA1**: SMR surface area, bin 1.
63. **SMR_VSA10**: SMR surface area, bin 10.
64. **SMR_VSA2**: SMR surface area, bin 2.
65. **SMR_VSA3**: SMR surface area, bin 3.
66. **SMR_VSA4**: SMR surface area, bin 4.
67. **SMR_VSA5**: SMR surface area, bin 5.
68. **SMR_VSA6**: SMR surface area, bin 6.
69. **SMR_VSA7**: SMR surface area, bin 7.
70. **SMR_VSA8**: SMR surface area, bin 8.
71. **SMR_VSA9**: SMR surface area, bin 9.
72. **SlogP_VSA1**: SlogP surface area, bin 1.
73. **SlogP_VSA10**: SlogP surface area, bin 10.
74. **SlogP_VSA11**: SlogP surface area, bin 11.
75. **SlogP_VSA12**: SlogP surface area, bin 12.
76. **SlogP_VSA2**: SlogP surface area, bin 2.
77. **SlogP_VSA3**: SlogP surface area, bin 3.
78. **SlogP_VSA4**: SlogP surface area, bin 4.
79. **SlogP_VSA5**: SlogP surface area, bin 5.
80. **SlogP_VSA6**: SlogP surface area, bin 6.
81. **SlogP_VSA7**: SlogP surface area, bin 7.
82. **SlogP_VSA8**: SlogP surface area, bin 8.
83. **SlogP_VSA9**: SlogP surface area, bin 9.
84. **TPSA**: Topological polar surface area, a predictor of drug absorption.
85. **EState_VSA1**: E-state surface area, bin 1.
86. **EState_VSA10**: E-state surface area, bin 10.
87. **EState_VSA11**: E-state surface area, bin 11.
88. **EState_VSA2**: E-state surface area, bin 2.
89. **EState_VSA3**: E-state surface area, bin 3.
90. **EState_VSA4**: E-state surface area, bin 4.
91. **EState_VSA5**: E-state surface area, bin 5.
92. **EState_VSA6**: E-state surface area, bin 6.
93. **EState_VSA7**: E-state surface area, bin 7.
94. **EState_VSA8**: E-state surface area, bin 8.
95. **EState_VSA9**: E-state surface area, bin 9.
96. **VSA_EState1**: VSA_E-state, bin 1.
97. **VSA_EState10**: VSA_E-state, bin 10.
98. **VSA_EState2**: VSA_E-state, bin 2.
99. **VSA_EState3**: VSA_E-state, bin 3.
100. **VSA_EState4**: VSA_E-state, bin 4.
101. **VSA_EState5**: VSA_E-state, bin 5.
102. **VSA_EState6**: VSA_E-state, bin 6.
103. **VSA_EState7**: VSA_E-state, bin 7.
104. **VSA_EState8**: VSA_E-state, bin 8.
105. **VSA_EState9**: VSA_E-state, bin 9.
106. **FractionCSP3**: Fraction of carbon atoms that are sp3 hybridized.
107. **HeavyAtomCount**: The number of heavy atoms (non-hydrogen atoms).
108. **NHOHCount**: The number of -OH and -NH groups.
109. **NOCount**: The number of nitrogen and oxygen atoms.
110. **NumAliphaticCarbocycles**: The number of aliphatic carbocycles.
111. **NumAliphaticHeterocycles**: The number of aliphatic heterocycles.
112. **NumAliphaticRings**: The number of aliphatic rings.
113. **NumAromaticCarbocycles**: The number of aromatic carbocycles.
114. **NumAromaticHeterocycles**: The number of aromatic heterocycles.
115. **NumAromaticRings**: The number of aromatic rings.
116. **NumHAcceptors**: The number of hydrogen bond acceptors.
117. **NumHDonors**: The number of hydrogen bond donors.
118. **NumHeteroatoms**: The number of heteroatoms (non-carbon and non-hydrogen atoms).
119. **NumRotatableBonds**: The number of rotatable bonds.
120. **NumSaturatedCarbocycles**: The number of saturated carbocycles.
121. **NumSaturatedHeterocycles**: The number of saturated heterocycles.
122. **NumSaturatedRings**: The number of saturated rings.
123. **RingCount**: The total number of rings in the molecule.
124. **MolLogP**: The logarithm of the partition coefficient (LogP) of the molecule.
125. **MolMR**: The molar refractivity of the molecule.

## Functional Group Descriptors
### Polar Functional Groups
1. **Hydroxyl Group ('[OH]')**: Indicates the presence of a hydroxyl group (-OH) in the molecule.
2. **Carbonyl Group ('C=O')**: Indicates the presence of a carbonyl group (C=O) in the molecule.
3. **Amide Group ('C(=O)N')**: Indicates the presence of an amide group (-C(=O)N-) in the molecule.
4. **Carboxyl Group ('C(=O)[OH]')**: Indicates the presence of a carboxyl group (-COOH) in the molecule.

### Non-Polar Functional Groups
5. **Alkyl ('[R]')**: Indicates the presence of an alkyl group (R-) in the molecule.
6. **Aromatic Rings ('c')**: Indicates the presence of aromatic rings in the molecule.
7. **Alkene ('C=C')**: Indicates the presence of an alkene group (C=C) in the molecule.

## Feature-Engineered Descriptors
1. **charge**: The total charge of the molecule.
2. **many_double_bonds**: The count of double bonds in the molecule.
3. **atoms_degree_0**: Number of atoms with zero degree (not connected to other atoms).
4. **atoms_degree_1**: Number of atoms with one connection.
5. **atoms_degree_2**: Number of atoms with two connections.
6. **atoms_degree_3**: Number of atoms with three connections.
7. **atoms_degree_4**: Number of atoms with four connections.
8. **atoms_degree_5**: Number of atoms with five connections.
9. **atoms_degree_6**: Number of atoms with six connections.
10. **atoms_valence_0**: Number of atoms with zero valence electrons.
11. **atoms_valence_1**: Number of atoms with one valence electron.
12. **atoms_valence_2**: Number of atoms with two valence electrons.
13. **atoms_valence_3**: Number of atoms with three valence electrons.
14. **atoms_valence_4**: Number of atoms with four valence electrons.
15. **atoms_valence_5**: Number of atoms with five valence electrons.
16. **atoms_valence_6**: Number of atoms with six valence electrons.
17. **atom_hybridization_S**: Number of atoms with S orbital hybridization.
18. **atom_hybridization_SP**: Number of atoms with SP orbital hybridization.
19. **atom_hybridization_SP2**: Number of atoms with SP2 orbital hybridization.
20. **atom_hybridization_SP3**: Number of atoms with SP3 orbital hybridization.
21. **atom_hybridization_SP3D**: Number of atoms with SP3D orbital hybridization.
22. **atom_hybridization_SP3D2**: Number of atoms with SP3D2 orbital hybridization.
23. **atom_hybridization_UNSPECIFIED**: Number of atoms with unspecified hybridization.
24. **aromatic_atoms**: Number of aromatic atoms in the molecule.
25. **single_bonds**: Number of single bonds in the molecule.
26. **double_bonds**: Number of double bonds in the molecule.
27. **triple_bonds**: Number of triple bonds in the molecule.
28. **aromatic_bonds**: Number of aromatic bonds in the molecule.
29. **zero_bonds**: Number of atoms with zero bonds.
30. **conjugated_bonds**: Number of conjugated bonds in the molecule.
31. **bonds_in_ring**: Number of bonds present in a ring structure.
32. **chirality_none**: Number of atoms without chirality.
33. **chirality_any**: Number of atoms with any chirality.
34. **chirality_z**: Number of atoms with Z chirality.
35. **chirality_e**: Number of atoms with E chirality.
36. **n_atoms**: Total number of atoms in the molecule.
37. **n_bonds**: Total number of bonds in the molecule.
38. **n_rings**: Total number of rings in the molecule.

## Descriptors Used by Sorkun
1. **nX**: Number of halogen atoms in the molecule.
2. **nHeavyAtom**: Number of non-hydrogen atoms in the molecule.
3. **nAromAtom**: Number of aromatic atoms in the molecule.
4. **nHBAcc**: Number of hydrogen bond acceptors.
5. **nHBDon**: Number of hydrogen bond donors.
6. **nRot**: Number of rotatable bonds.
7. **nBonds**: Total number of bonds in the molecule.
8. **nAromBond**: Number of aromatic bonds in the molecule.
9. **nBondsO**: Number of bonds to oxygen atoms.
10. **nBondsS**: Number of bonds to sulfur atoms.
11. **nBondsD**: Number of double bonds in the molecule.
12. **nBondsT**: Number of triple bonds in the molecule.
13. **VMcGowan**: McGowan volume, a descriptor related to molecular size.
14. **TopoPSA(NO)**: Topological polar surface area, calculated considering only oxygen and nitrogen atoms.
15. **TopoPSA**: Topological polar surface area, considering all polar atoms.
16. **LabuteASA**: Labute’s Approximation to Molecular Surface Area.
17. **apol**: Sum of atomic polarizabilities (including implicit hydrogens).
18. **bpol**: Sum of atomic polarizabilities (excluding hydrogens).
19. **nAcid**: Number of acidic groups in the molecule.
20. **nBase**: Number of basic groups in the molecule.
21. **ECIndex**: Eccentric connectivity index, a topological descriptor.
22. **GGI1**: Gutman Molecular Topological Index.
23. **JGI1**: Topological charge index of the first order.
24. **SLogP**: Logarithm of the partition coefficient (SLogP).
25. **SMR**: Sum of atomic molar refractivity values.
26. **BertzCT**: Bertz complexity index, related to the molecular complexity.
27. **BalabanJ**: Balaban’s J index, a topological descriptor.
28. **WPol**: Weighted path order 3.
29. **Zagreb1**: First Zagreb index, a topological descriptor.
30. **ABC**: Aromatic bond count.
31. **ABCGG**: Generalized graph approach to aromatic bond count.
32. **nRing**: Number of rings in the molecule.
33. **nHRing**: Number of rings with hydrogen atoms.
34. **naRing**: Number of aromatic rings in the molecule.
35. **naHRing**: Number of aromatic rings with hydrogen atoms.
36. **nARing**: Number of rings with atoms other than carbon.
37. **nFRing**: Number of fused rings in the molecule.
38. **NsCH3**: Number of sulfur-carbon single bonds where carbon is connected to three other atoms.
39. **NdCH2**: Number of carbon atoms doubly bonded to another atom and singly bonded to two other atoms.
40. **NssCH2**: Number of carbon atoms singly bonded to a sulfur atom and two other atoms.
41. **NtCH**: Number of carbon atoms triply bonded to another atom and singly bonded to one other atom.
42. **NdsCH**: Number of carbon atoms doubly bonded to a sulfur atom and singly bonded to one other atom.
43. **NaaCH**: Number of carbon atoms singly bonded to two other atoms and a sulfur atom.
44. **NsssCH**: Number of carbon atoms singly bonded to three other atoms and a sulfur atom.
45. **NddC**: Number of carbon atoms doubly bonded to two other atoms.
46. **NtsC**: Number of carbon atoms triply bonded to another atom and singly bonded to one other atom.
47. **NdssC**: Number of carbon atoms doubly bonded to a sulfur atom and singly bonded to one other atom.
48. **NaasC**: Number of carbon atoms singly bonded to a sulfur atom and one other atom.
49. **NaaaC**: Number of carbon atoms singly bonded to three other atoms.
50. **NssssC**: Number of carbon atoms singly bonded to four other atoms.
51. **NsNH2**: Number of nitrogen atoms singly bonded to one other atom and two hydrogens.
52. **NssNH**: Number of nitrogen atoms singly bonded to two other atoms and one hydrogen.
53. **NaaNH**: Number of nitrogen atoms singly bonded to two other atoms and one hydrogen.
54. **NtN**: Number of nitrogen atoms triply bonded to another atom.
55. **NdsN**: Number of nitrogen atoms doubly bonded to another atom and singly bonded to one other atom.
56. **NaaN**: Number of nitrogen atoms singly bonded to two other atoms.
57. **NsssN**: Number of nitrogen atoms singly bonded to three other atoms.
58. **NaasN**: Number of nitrogen atoms singly bonded to two other atoms.
59. **NsOH**: Number of oxygen atoms singly bonded to one other atom and one hydrogen.
60. **NdO**: Number of oxygen atoms doubly bonded to another atom.
61. **NssO**: Number of oxygen atoms singly bonded to two other atoms.
62. **NaaO**: Number of oxygen atoms singly bonded to two other atoms.
63. **NsF**: Number of fluorine atoms singly bonded to one other atom.
64. **NdsssP**: Number of phosphorus atoms doubly bonded to another atom and singly bonded to three other atoms.
65. **NdS**: Number of sulfur atoms doubly bonded to another atom.
66. **NssS**: Number of sulfur atoms singly bonded to two other atoms.
67. **NaaS**: Number of sulfur atoms singly bonded to two other atoms.
68. **NdssS**: Number of sulfur atoms doubly bonded to another atom and singly bonded to one other atom.
69. **NddssS**: Number of sulfur atoms doubly bonded to two other atoms.
70. **NsCl**: Number of chlorine atoms singly bonded to one other atom.
71. **NsBr**: Number of bromine atoms singly bonded to one other atom.
72. **NsI**: Number of iodine atoms singly bonded to one other atom.
73. **SsCH3**: Number of carbon atoms singly bonded to one sulfur atom and three other atoms.
74. **SdCH2**: Number of carbon atoms doubly bonded to another atom and singly bonded to two other atoms.
75. **SssCH2**: Number of carbon atoms singly bonded to two other atoms and one sulfur atom.
76. **StCH**: Number of carbon atoms triply bonded to another atom and singly bonded to one other atom.
77. **SdsCH**: Number of carbon atoms doubly bonded to a sulfur atom and singly bonded to one other atom.
78. **SaaCH**: Number of carbon atoms singly bonded to a sulfur atom and two other atoms.
79. **SsssCH**: Number of carbon atoms singly bonded to three other atoms and a sulfur atom.
80. **SddC**: Number of carbon atoms doubly bonded to two other atoms.
81. **StsC**: Number of carbon atoms triply bonded to another atom and singly bonded to one other atom.
82. **SdssC**: Number of carbon atoms doubly bonded to a sulfur atom and singly bonded to one other atom.
83. **SaasC**: Number of carbon atoms singly bonded to two other atoms and a sulfur atom.
84. **SaaaC**: Number of carbon atoms singly bonded to three other atoms.
85. **SssssC**: Number of carbon atoms singly bonded to four other atoms.
86. **SsNH2**: Number of nitrogen atoms singly bonded to one other atom and two hydrogens.
87. **SssNH**: Number of nitrogen atoms singly bonded to two other atoms and one hydrogen.
88. **SaaNH**: Number of nitrogen atoms singly bonded to two other atoms and one hydrogen.
89. **StN**: Number of nitrogen atoms triply bonded to another atom.
90. **SdsN**: Number of nitrogen atoms doubly bonded to another atom and singly bonded to one other atom.
91. **SaaN**: Number of nitrogen atoms singly bonded to two other atoms.
92. **SsssN**: Number of nitrogen atoms singly bonded to three other atoms.
93. **SaasN**: Number of nitrogen atoms singly bonded to two other atoms.
94. **SsOH**: Number of oxygen atoms singly bonded to one other atom and one hydrogen.
95. **SdO**: Number of oxygen atoms doubly bonded to another atom.
96. **SssO**: Number of oxygen atoms singly bonded to two other atoms.
97. **SaaaO**: Number of oxygen atoms singly bonded to two other atoms.
98. **SsF**: Number of fluorine atoms singly bonded to one other atom.
99. **SdsssP**: Number of phosphorus atoms doubly bonded to another atom and singly bonded to three other atoms.
100. **SdS**: Number of sulfur atoms doubly bonded to another atom.
101. **SssS**: Number of sulfur atoms singly bonded to two other atoms.
102. **SaaaS**: Number of sulfur atoms singly bonded to two other atoms.
103. **SdssS**: Number of sulfur atoms doubly bonded to another atom and singly bonded to one other atom.
104. **SddssS**: Number of sulfur atoms doubly bonded to two other atoms.
105. **SsCl**: Number of chlorine atoms singly bonded to one other atom.
106. **SsBr**: Number of bromine atoms singly bonded to one other atom.
107. **SsI**: Number of iodine atoms singly bonded to one other atom.
108. **C**: Total number of carbon atoms.
109. **H**: Total number of hydrogen atoms.
110. **Br**: Total number of bromine atoms.
111. **N**: Total number of nitrogen atoms.
112. **O**: Total number of oxygen atoms.
113. **I**: Total number of iodine atoms.
114. **Cl**: Total number of chlorine atoms.
115. **S**: Total number of sulfur atoms.
116. **F**: Total number of fluorine atoms.
117. **P**: Total number of phosphorus atoms.
118. **As**: Total number of arsenic atoms.
119. **Si**: Total number of silicon atoms.
120. **Se**: Total number of selenium atoms.
121. **Sn**: Total number of tin atoms.
122. **Hg**: Total number of mercury atoms.
123. **Ge**: Total number of germanium atoms.
