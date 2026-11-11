"""
Feature & analysis descriptions for PROFETdb right-click info popups.

You only need to fill in FEATURE_BRIEFS and ANALYSIS_BRIEFS below.
URLs and wiring are handled for you.
"""

FEATURE_PAGE = "https://tools.tartaglialab.com/profetDB/featuredescription"
DOCUMENTATION_PAGE = "https://tools.tartaglialab.com/profetDB/documentation"

def feature_url(feature_id: str) -> str:
    # If the page expects a query param instead of an anchor, switch to:
    #   return f"{FEATURE_PAGE}?feature={feature_id}"
    #return f"{FEATURE_PAGE}#{feature_id}" # do this if there are anchors in page
    return f"{FEATURE_PAGE}"

def analysis_url(slug: str) -> str:
    #return f"{DOCUMENTATION_PAGE}#{slug}" # do this if there are anchors in page
    return f"{DOCUMENTATION_PAGE}"
    
PHYSICAL_CHEMICAL = ['BetaSheet_1', 'Aggregation_1', 'Disorder_1', 'Aggregation_2', 'Membrane_1', 'Hydrophobicity_1', 'BetaSheet_2', 'BetaSheet_3', 'AlphaHelix_1', 'Membrane_2', 'BetaSheet_4', 'Burial_1', 'Burial_2', 'Hydrophobicity_2', 'NucleicAcidBinding_1', 'AlphaHelix_2', 'BetaSheet_5', 'charge', 'LowFolding_1', 'Membrane_3', 'Turn_1', 'fg', 'BetaSheet_6', 'NucleicAcidBinding_2', 'BetaSheet_7', 'Disorder_2', 'Membrane_4', 'Burial_3', 'Hydrophobicity_3', 'BetaSheet_8', 'Burial_4', 'Disorder_3', 'Aggregation_3', 'Membrane_5', 'Hydrophobicity_4', 'Hydrophobicity_5', 'Hydrophobicity_6', 'Aggregation_4', 'Aggregation_5', 'Disorder_4', 'Aggregation_6', 'Hydrophobicity_7', 'Membrane_6', 'AlphaHelix_3', 'Disorder_5', 'Disorder_6', 'AlphaHelix_4', 'Burial_5', 'NucleicAcidBinding_3', 'AlphaHelix_5', 'NucleicAcidBinding_4', 'BetaSheet_9', 'Aggregation_7', 'Hydrophobicity_8', 'rg', 'Disorder_7', 'NucleicAcidBinding_5', 'Disorder_8', 'Disorder_9', 'Aggregation_8', 'Burial_6', 'AlphaHelix_6', 'NucleicAcidBinding_6', 'AlphaHelix_7', 'Aggregation_9', 'BetaSheet_10', 'Burial_7', 'Disorder_10', 'Hydrophobicity_9', 'AlphaHelix_8', 'Burial_8', 'AlphaHelix_9', 'NucleicAcidBinding_7', 'NucleicAcidBinding_8', 'NucleicAcidBinding_9', 'AlphaHelix_10', 'NucleicAcidBinding_10', 'Hydrophobicity_10', 'Membrane_7', 'Membrane_8', 'Membrane_9', 'Burial_9']
STRUCTURAL = ['n_contacts', 'RG_protein', 'rmsd_prot', 'Length', 'fullCharge', 'extCharge', 'Coil_FullSeq', 'AlphaHelix_FullSeq', 'BetaBridge_FullSeq', 'Strand_FullSeq', '3_10Helix_FullSeq', 'PiHelix_FullSeq', 'Turn_FullSeq', 'Bend_FullSeq', 'Coil_ExtSeq', 'AlphaHelix_ExtSeq', 'BetaBridge_ExtSeq', 'Strand_ExtSeq', '3_10Helix_ExtSeq', 'PiHelix_ExtSeq', 'Turn_ExtSeq', 'Bend_ExtSeq', 'asa_mean', 'asa_std', 'RBD_Full_min_0', 'RBD_Full_min_1', 'RBD_Full_min_2', 'RBD_int_min_0', 'RBD_int_min_1', 'RBD_int_min_2', 'RBD_ext_min_0', 'RBD_ext_min_1', 'RBD_ext_min_2', 'average_plddt', 'stddev_plddt', 'n_contacts_norm', 'RG_protein_norm']
GRAPH = ['number_of_nodes', 'number_of_edges', 'Density', 'size_of_largest_connected_component', 'avgerage_degree', 'degree_assortativity', 'betweeness_centrality', 'eigenvector_centrality', 'closeness_centrality', 'katz_centrality', 'harmonic_centrality', 'pagerank', 'hub_scores', 'authority_scores', 'local_clustering_coefficient', 'average_clustering_coefficient', 'global_transitivity', 'degree_centrality', 'average_shortest_path_length', 'diameter', 'entropy']
THERMOSTABILITY = ['Ts_tot_PN', 'Ts_tot_N', 'Ts_tot_P', 'Ts_C_PN', 'Ts_C_N', 'Ts_C_P', 'Ts_lj_PN', 'Ts_lj_N', 'Ts_lj_P']
EXPRESSION = ['mean_expression','tau']
PROFILE_ONLY = ['Membrane_8', 'Aggregation_4', 'Aggregation_3', 'LowFolding_1', 'Aggregation_5', 'Aggregation_2', 'Aggregation_7', 'Aggregation_6', 'Aggregation_8', 'Aggregation_9', 'Aggregation_1', 'AlphaHelix_8', 'Turn_1', 'Disorder_8', 'BetaSheet_2', 'AlphaHelix_9', 'BetaSheet_10', 'Burial_1', 'Disorder_5', 'BetaSheet_8', 'charge', 'Disorder_9', 'Membrane_2', 'Disorder_10', 'Burial_4', 'fg', 'Hydrophobicity_8', 'Hydrophobicity_3', 'Hydrophobicity_10', 'Hydrophobicity_5', 'Hydrophobicity_6', 'Hydrophobicity_7', 'Hydrophobicity_4', 'Hydrophobicity_1', 'Hydrophobicity_9', 'Hydrophobicity_2', 'AlphaHelix_5', 'BetaSheet_6', 'Membrane_9', 'Membrane_4', 'Burial_2', 'Burial_8', 'Membrane_1', 'Membrane_5', 'AlphaHelix_3', 'AlphaHelix_7', 'AlphaHelix_6', 'rg', 'AlphaHelix_10', 'AlphaHelix_4', 'BetaSheet_3', 'BetaSheet_4', 'BetaSheet_7', 'BetaSheet_5', 'BetaSheet_1', 'Disorder_1', 'Disorder_2', 'AlphaHelix_1', 'Disorder_7', 'NucleicAcidBinding_2', 'NucleicAcidBinding_9', 'NucleicAcidBinding_8', 'NucleicAcidBinding_4', 'NucleicAcidBinding_1', 'NucleicAcidBinding_5', 'NucleicAcidBinding_7', 'NucleicAcidBinding_6', 'NucleicAcidBinding_10', 'NucleicAcidBinding_3', 'Burial_5', 'Burial_7', 'Burial_9', 'Burial_6', 'Burial_3', 'AlphaHelix_2', 'BetaSheet_9', 'Disorder_6', 'Disorder_3', 'Membrane_7', 'Membrane_3', 'Membrane_6', 'Disorder_4', 'local_clustering_coefficient', 'degrees', 'degree_centrality', 'betweeness_centrality', 'eigenvector_centrality', 'closeness_centrality', 'katz_centrality', 'harmonic_centrality', 'pagerank', 'hub_scores', 'authority_scores', 'eccentricity', 'triangles', 'eigenvalues']

ALL_FEATURE_IDS = (PHYSICAL_CHEMICAL + STRUCTURAL + GRAPH + THERMOSTABILITY + EXPRESSION + PROFILE_ONLY)

# CHANGE THIS IN LIVE
FEATURE_BRIEFS = {
    'Membrane_1': 'A ranking of the amino acids by their preference to lie in lipid contact (membrane-buried), obtained from a structural-prediction algorithm for membrane-bound proteins calibrated against bacteriorhodopsin.',
    'Membrane_2': "Amino-acid composition (percentage) characteristic of integral-membrane proteins, from Cedano et al.'s correlation analysis linking residue composition to protein cellular location.",
    'Membrane_3': 'Amino-acid composition of the transmembrane regions of non-mitochondrial membrane proteins, from a compositional analysis showing membrane proteins have a distinctly hydrophobic residue content.',
    'Membrane_4': 'Knowledge-based membrane-propensity scale derived from transmembrane-helix segments in membrane-protein crystal structures (3D topology databases), optimized so membrane segments are energetically favoured over typical soluble sequences.',
    'Membrane_5': 'Normalized amino-acid composition of membrane proteins, capturing the hydrophobic-residue preferences that distinguish membrane from soluble proteins.',
    'Membrane_6': "Transmembrane-helix turn propensity measuring each residue's ability to form a tight turn within a hydrophobic membrane segment; values correlate closely with residue hydrophobicity.",
    'Membrane_7': 'Amino-acid composition of the transmembrane regions of mitochondrial membrane proteins, which are notably rich in hydrophobic residues.',
    'Membrane_8': 'Overall amino-acid composition of membrane proteins, used to discriminate membrane from soluble sequences by their hydrophobic-residue content.',
    'Membrane_9': 'Knowledge-based membrane-propensity scale derived from 1D helix-in-membrane-protein topology databases, an alternative parameterization to the 3D-structure-based scale.',
    'Turn_1': 'Averaged transmembrane-helix turn propensity per residue, from experimental measurement of turn-forming ability within a model hydrophobic segment.',
    'NucleicAcidBinding_1': 'Propensity of a residue to sit at the centre of an RNA-binding interface, derived from the RNABindR analysis of experimentally determined protein-RNA complexes.',
    'NucleicAcidBinding_2': 'Enrichment of the residue among proteins bearing classical RNA-binding domains, from the mRNA-interactome atlas of mammalian RNA-binding proteins.',
    'NucleicAcidBinding_3': 'Propensity of a residue to lie immediately upstream of an RNA-binding interface residue, from the RNABindR interface analysis.',
    'NucleicAcidBinding_4': 'Residue propensity to occur at protein-RNA interfaces, tabulated from the PRIDB database of protein-RNA interfaces extracted from PDB complexes.',
    'NucleicAcidBinding_5': 'Propensity of a residue to lie immediately downstream of an RNA-binding interface residue, from the RNABindR interface analysis.',
    'NucleicAcidBinding_6': 'Enrichment of the residue among proteins bearing non-classical RNA-binding domains, from the mammalian mRNA-interactome atlas.',
    'NucleicAcidBinding_7': 'Enrichment of the residue among the experimentally captured mRNA interactome (all identified RNA-binding proteins).',
    'NucleicAcidBinding_8': 'Residue propensity at protein-RNA interfaces mediated by water (HOH)-bridged contacts, from the PRIDB interface database.',
    'NucleicAcidBinding_9': 'Residue propensity at protein-RNA interfaces forming direct hydrogen-bond (HB) contacts, from the PRIDB interface database.',
    'NucleicAcidBinding_10': 'Enrichment of the residue among RNA-binding proteins that lack a recognized/known RNA-binding domain, from the mammalian mRNA-interactome atlas.',
    'Disorder_1': 'Normalized frequency with which each residue occurs in random-coil conformation, from a statistical analysis of 95 proteins for secondary-structure prediction.',
    'Disorder_2': 'Normalized frequency of the coil (irregular/chain-reversal) conformation per residue, derived from a statistical treatment of protein-folding conformations.',
    'Disorder_3': 'TOP-IDP intrinsic-disorder propensity, an amino-acid scale optimized to best discriminate ordered from intrinsically disordered residues.',
    'Disorder_4': 'Residue order/disorder propensity based on an unfolded-over-folded preference, one of the disorder-discriminating attributes evaluated in the TOP-IDP study.',
    'Disorder_5': 'Residue disorder propensity based on crystallographic B-factor (flexibility) values, among the scales assessed in the TOP-IDP analysis.',
    'Disorder_6': 'The Chou-Fasman conformational parameter for the coil state, related to residue steric and intermolecular-force properties.',
    'Disorder_7': 'Normalized relative frequency of the coil conformation per residue, from a statistical characterization of bends and chain reversals in proteins.',
    'Disorder_8': 'Average residue flexibility index quantifying the positional (fluctuational) mobility of each amino acid in globular proteins.',
    'Disorder_9': 'Coil-conformation preference parameter from the class-based double-prediction secondary-structure method.',
    'Disorder_10': 'Disorder propensity derived from residue frequencies in experimentally verified disordered regions of the DisProt database, as used in the TOP-IDP analysis.',
    'AlphaHelix_1': 'Normalized frequency of each residue in the middle positions of alpha-helices, from a conformational analysis of globular proteins (the study that also defined the reverse turn).',
    'AlphaHelix_2': 'Relative frequency of the residue in alpha-helical conformation, from a compilation of conformational-property distributions in signal and nascent peptides.',
    'AlphaHelix_3': 'Normalized alpha-helix frequency parameter from a four-state (helix/extended/bend/coil) conformational analysis of proteins.',
    'AlphaHelix_4': 'Normalized alpha-helix propensity computed specifically from proteins of the all-alpha structural class, using class-specific preference parameters.',
    'AlphaHelix_5': 'Information-theoretic alpha-helix preference measure per residue, from a statistical sequence-conformation analysis.',
    'AlphaHelix_6': 'Normalized (unweighted) frequency of the residue in alpha-helix, from an analysis of secondary-structure occurrence across more than 50 globular proteins.',
    'AlphaHelix_7': "The Chou-Fasman P-alpha conformational parameter: each residue's normalized frequency of occurrence in alpha-helix, tabulated from X-ray protein structures for their secondary-structure prediction method.",
    'AlphaHelix_8': 'Alpha-helix preference parameter from the class-based double-prediction secondary-structure method.',
    'AlphaHelix_9': 'Average relative probability of alpha-helix as a function of local hydrophobicity, from an analysis linking local hydrophobicity to secondary-structure formation.',
    'AlphaHelix_10': 'Weighted normalized frequency of the residue in alpha-helix, from a survey of secondary structure across globular proteins.',
    'Aggregation_1': 'Intrinsic aggregation propensity (medium-sensitivity setting) predicting amyloid-prone regions while accounting for competition with folding, for structured proteins.',
    'Aggregation_2': 'Beta-aggregation propensity (high setting) from the TANGO algorithm, based on the physico-chemical principles of beta-sheet formation with fully buried aggregate cores.',
    'Aggregation_3': 'Intrinsic per-residue aggregation propensity used to identify aggregation-prone and aggregation-susceptible regions in neurodegenerative disease proteins.',
    'Aggregation_4': 'Aggregation propensity from AGGRESCAN, based on an in-vivo-derived per-residue aggregation-propensity scale for predicting aggregation-prone segments.',
    'Aggregation_5': 'Beta-aggregation propensity from a physicochemical model predicting aggregation rate and aggregation-prone segments in polypeptide sequences.',
    'Aggregation_6': 'Beta-aggregation propensity (low setting) from the TANGO algorithm based on beta-sheet formation principles.',
    'Aggregation_7': 'Intrinsic aggregation propensity (high-sensitivity setting) predicting amyloid-prone regions while accounting for folding competition.',
    'Aggregation_8': 'Intrinsic aggregation propensity (low-sensitivity setting) predicting amyloid-prone regions while accounting for folding competition.',
    'Aggregation_9': 'Beta-aggregation propensity (medium setting) from the TANGO algorithm based on beta-sheet formation principles.',
    'LowFolding_1': 'Propensity for low folding (equivalently high aggregation), capturing the interplay between folding and aggregation via the CamFold/Zyggregator approach.',
    'BetaSheet_1': 'Normalized beta-sheet propensity computed from proteins of the all-beta structural class, using class-specific preference parameters.',
    'BetaSheet_2': 'Average relative probability of beta-sheet as a function of local hydrophobicity, from an analysis linking local hydrophobicity to secondary structure.',
    'BetaSheet_3': "The Chou-Fasman P-beta conformational parameter: each residue's normalized frequency of occurrence in beta-sheet, tabulated from X-ray protein structures for their secondary-structure prediction method.",
    'BetaSheet_4': 'Normalized frequency of each residue in beta-sheet, from a conformational analysis of globular proteins (the study that also defined the reverse turn).',
    'BetaSheet_5': 'Weighted normalized frequency of the residue in beta-sheet, from a survey of secondary structure across globular proteins.',
    'BetaSheet_6': 'Information-theoretic beta-(pleated-)sheet preference measure per residue, from a statistical sequence-conformation analysis.',
    'BetaSheet_7': 'Unweighted normalized frequency of the residue in beta-sheet, from a survey of secondary structure across globular proteins.',
    'BetaSheet_8': 'Beta-sheet preference parameter from the class-based double-prediction secondary-structure method.',
    'BetaSheet_9': 'Relative frequency of the residue in beta-sheet conformation, from a compilation of conformational-property distributions in signal and nascent peptides.',
    'BetaSheet_10': 'Average relative probability of the residue in inner beta-sheet strands as a function of local hydrophobicity.',
    'Burial_1': 'Average volume occupied by each residue when buried in the protein interior, reflecting the close, crystal-like packing of protein cores.',
    'Burial_2': 'Mean solvent-accessible surface area each residue buries upon folding, correlated with residue hydrophobicity.',
    'Burial_3': 'Ratio of buried to accessible molar fractions for each residue, from a solvent-accessibility analysis distinguishing interior from surface residues.',
    'Burial_4': "Free energy of transferring each residue from the aqueous exterior to the buried protein interior (95%-buried reference), based on water-to-cyclohexane side-chain distribution coefficients that closely track residues' inside/outside distribution in globular proteins.",
    'Burial_5': 'Percentage of each residue type found buried (solvent-inaccessible) in globular proteins, from an analysis of side-chain conformations and accessibility across 19 protein structures.',
    'Burial_6': 'Proportion of each residue type that is at least 95% buried in the protein interior, from an analysis of accessible and buried surface areas in proteins.',
    'Burial_7': 'Propensity of each residue to be located inside rather than on the surface of a protein, from an inside/outside preference analysis.',
    'Burial_8': 'Mean volume occupied by each residue buried in protein interiors, reflecting the exceptionally high packing density of the protein core.',
    'Burial_9': 'Proportion of each residue type that is essentially fully (100%) buried in the protein interior, from an analysis of accessible and buried surface areas in proteins.',
    'Hydrophobicity_1': 'Membrane-buried conformational-preference (hydrophobicity) parameter for predicting transmembrane helices in integral membrane proteins.',
    'Hydrophobicity_2': 'Optimal-matching hydrophobicity (OMH) value per residue, derived from observed amino-acid replacement frequencies among proteins of similar three-dimensional structure so that sequence hydrophobicity patterns track structural similarity.',
    'Hydrophobicity_3': 'Hydrophobicity parameters covering the standard residues (and post-/co-translationally modified ones), validated on hydrophobicity profiles and moments.',
    'Hydrophobicity_4': 'The Kyte-Doolittle hydropathy value per residue, widely used with a sliding window to locate buried cores and membrane-spanning segments.',
    'Hydrophobicity_5': 'Normalized consensus hydrophobicity scale, used with the hydrophobic-moment plot to characterize membrane and surface helices.',
    'Hydrophobicity_6': 'Hansch hydrophobicity parameters (pi) for the amino-acid side chains, measured from octanol-water partitioning of their N-acetyl-amino-acid amides; the values correlate tightly with the degree of side-chain solvent exposure in globular proteins.',
    'Hydrophobicity_7': 'Residue hydrophobicity expressed as the buried/accessible preference from a solvent-accessibility analysis of globular proteins.',
    'Hydrophobicity_8': 'Side-chain hydrophobicity (partition coefficient) computed with an extended fragment method for the amino acids.',
    'Hydrophobicity_9': 'Hydrophilicity/hydrophobicity scale based on water/organic partitioning of N-acetyl amino-acid amides, accounting for the reduction of side-chain hydrophilicity by flanking peptide bonds.',
    'Hydrophobicity_10': 'Hydrophobicity scale derived from the free energy of transfer of each amino-acid residue from solution to the air-water surface, measured via the surface tension of amino-acid solutions.',
    'rg': 'Local arginine-glycine (RG) composition of the sequence, computed per residue and as a whole-protein mean.',
    'fg': 'Local phenylalanine-glycine (FG) composition of the sequence, computed per residue and as a whole-protein mean.',
    'charge': 'Net electrostatic charge of the protein derived from its amino-acid composition.',
    'n_contacts': 'Number of residue-residue contacts in the AlphaFold 3D structure, counting amino-acid pairs closer than 5 Angstrom.',
    'RG_protein': 'Radius of gyration of the protein: the root-mean-square distance of its atoms from the centre of mass, describing overall compactness.',
    'rmsd_prot': "Average pairwise distance between residues in the 3D structure, a measure of the protein's spatial spread.",
    'Length': 'Total number of residues in the protein sequence.',
    'fullCharge': 'Net electrostatic charge of the whole protein, summed over all residues.',
    'extCharge': 'Net electrostatic charge of the solvent-exposed residues only, computed from the 3D structure.',
    'Coil_FullSeq': 'Percentage of residues not assigned to any regular secondary-structure element (coil/loop) over the full structure; used as a proxy for flexible/disordered content.',
    'AlphaHelix_FullSeq': 'Percentage of residues assigned to alpha-helix (DSSP class H) over the full 3D structure.',
    'BetaBridge_FullSeq': 'Percentage of residues assigned to isolated beta-bridge (DSSP class B) over the full 3D structure.',
    'Strand_FullSeq': 'Percentage of residues assigned to extended beta-strand (DSSP class E) over the full 3D structure.',
    '3_10Helix_FullSeq': 'Percentage of residues assigned to 3-10 helix (DSSP class G) over the full 3D structure.',
    'PiHelix_FullSeq': 'Percentage of residues assigned to pi-helix (DSSP class I) over the full 3D structure.',
    'Turn_FullSeq': 'Percentage of residues assigned to hydrogen-bonded turn (DSSP class T) over the full 3D structure.',
    'Bend_FullSeq': 'Percentage of residues assigned to bend (DSSP class S) over the full 3D structure.',
    'Coil_ExtSeq': 'Percentage of solvent-exposed residues not assigned to any regular secondary-structure element (coil/loop).',
    'AlphaHelix_ExtSeq': 'Percentage of solvent-exposed residues assigned to alpha-helix (DSSP class H).',
    'BetaBridge_ExtSeq': 'Percentage of solvent-exposed residues assigned to isolated beta-bridge (DSSP class B).',
    'Strand_ExtSeq': 'Percentage of solvent-exposed residues assigned to extended beta-strand (DSSP class E).',
    '3_10Helix_ExtSeq': 'Percentage of solvent-exposed residues assigned to 3-10 helix (DSSP class G).',
    'PiHelix_ExtSeq': 'Percentage of solvent-exposed residues assigned to pi-helix (DSSP class I).',
    'Turn_ExtSeq': 'Percentage of solvent-exposed residues assigned to hydrogen-bonded turn (DSSP class T).',
    'Bend_ExtSeq': 'Percentage of solvent-exposed residues assigned to bend (DSSP class S).',
    'asa_mean': 'Mean solvent-accessible surface area (ASA) over all residues, computed from the 3D structure.',
    'asa_std': 'Standard deviation of the per-residue solvent-accessible surface area, capturing heterogeneity of residue exposure.',
    'RBD_Full_min_0': 'Highest structural-similarity score between the protein and experimentally resolved RNA-binding proteins/domains, over the full structure.',
    'RBD_Full_min_1': 'Second-highest structural-similarity score to experimentally resolved RNA-binding proteins, over the full structure.',
    'RBD_Full_min_2': 'Third-highest structural-similarity score to experimentally resolved RNA-binding proteins, over the full structure.',
    'RBD_int_min_0': 'Highest structural-similarity score to experimentally resolved RNA-binding proteins, computed over the internal (buried) sequence.',
    'RBD_int_min_1': 'Second-highest structural-similarity score to experimentally resolved RNA-binding proteins, over the internal (buried) sequence.',
    'RBD_int_min_2': 'Third-highest structural-similarity score to experimentally resolved RNA-binding proteins, over the internal (buried) sequence.',
    'RBD_ext_min_0': 'Highest structural-similarity score to experimentally resolved RNA-binding proteins, computed over the solvent-exposed sequence.',
    'RBD_ext_min_1': 'Second-highest structural-similarity score to experimentally resolved RNA-binding proteins, over the solvent-exposed sequence.',
    'RBD_ext_min_2': 'Third-highest structural-similarity score to experimentally resolved RNA-binding proteins, over the solvent-exposed sequence.',
    'average_plddt': 'Mean per-residue pLDDT confidence score of the AlphaFold model; consistently low pLDDT regions tend to correspond to intrinsic disorder.',
    'stddev_plddt': 'Standard deviation of the per-residue pLDDT confidence scores across the structure.',
    'n_contacts_norm': 'Number of residue contacts (pairs closer than 5 Angstrom) normalized by protein length, giving a size-independent contact density.',
    'RG_protein_norm': 'Radius of gyration normalized by protein length, giving a size-independent measure of compactness.',
    'number_of_nodes': 'Number of nodes in the Residue Interaction Network (RIN), equal to the number of residues in the structure.',
    'number_of_edges': 'Number of edges in the RIN, i.e. the count of residue-residue contacts closer than 7 Angstrom.',
    'Density': 'Density of the RIN: the ratio of present edges to the maximum possible number of edges among its nodes.',
    'size_of_largest_connected_component': 'Number of residues in the largest connected component of the RIN.',
    'avgerage_degree': 'Average node degree of the RIN, i.e. the mean number of contacts per residue.',
    'degree_assortativity': 'Degree assortativity of the RIN: the Pearson correlation between the degrees of connected residues, indicating whether high-degree residues preferentially contact other high-degree residues.',
    'degree_centrality': 'Degree centrality per residue: the node degree normalized by the maximum possible degree in the RIN.',
    'betweeness_centrality': 'Betweenness centrality per residue: the fraction of all-pairs shortest paths in the RIN that pass through the residue, highlighting bottleneck residues.',
    'eigenvector_centrality': "Eigenvector centrality per residue: influence measured so that a residue's centrality is proportional to the sum of its neighbours' centralities in the RIN.",
    'closeness_centrality': 'Closeness centrality per residue: the reciprocal of the mean shortest-path distance from the residue to all others in the RIN.',
    'katz_centrality': 'Katz centrality per residue: centrality summing over all paths reaching a residue in the RIN, with longer paths attenuated by a decay factor.',
    'harmonic_centrality': 'Harmonic centrality per residue: the sum of the reciprocals of shortest-path distances from the residue to all others in the RIN.',
    'pagerank': 'PageRank per residue: the stationary probability of a random walk on the RIN, ranking residues by connectivity importance.',
    'hub_scores': 'HITS hub score per residue in the RIN: high for residues connected to many high-authority residues.',
    'authority_scores': 'HITS authority score per residue in the RIN: high for residues connected to many high-hub residues.',
    'local_clustering_coefficient': "Local clustering coefficient per residue: the fraction of a residue's RIN neighbours that are themselves in contact, measuring local cohesion.",
    'average_clustering_coefficient': 'Mean of the local clustering coefficients over all residues in the RIN.',
    'global_transitivity': 'Global transitivity of the RIN: three times the number of triangles divided by the number of connected residue triples.',
    'average_shortest_path_length': 'Average shortest-path length of the RIN: the mean number of contact steps between all residue pairs.',
    'diameter': 'Diameter of the RIN: the longest shortest-path distance between any two residues.',
    'entropy': "Shannon entropy of the RIN's vertex-degree distribution, quantifying the heterogeneity of residue connectivity.",
    'Ts_tot_PN': 'Total interaction energy (combined Coulombic + Lennard-Jones) summed over polar and non-polar residue interactions in the structure, from the energy-weighted interaction-graph representation used as a proxy for thermal stability.',
    'Ts_tot_N': 'Total interaction energy (combined Coulombic + Lennard-Jones) over non-polar residue interactions only.',
    'Ts_tot_P': 'Total interaction energy (combined Coulombic + Lennard-Jones) over polar residue interactions only.',
    'Ts_C_PN': 'Coulombic (electrostatic) interaction energy summed over polar and non-polar residue interactions in the structure.',
    'Ts_C_N': 'Coulombic (electrostatic) interaction energy over non-polar residue interactions only.',
    'Ts_C_P': 'Coulombic (electrostatic) interaction energy over polar residue interactions only.',
    'Ts_lj_PN': 'Lennard-Jones (van der Waals) interaction energy summed over polar and non-polar residue interactions in the structure.',
    'Ts_lj_N': 'Lennard-Jones (van der Waals) interaction energy over non-polar residue interactions only.',
    'Ts_lj_P': 'Lennard-Jones (van der Waals) interaction energy over polar residue interactions only.',
    'mean_expression': 'Mean protein abundance across all available human tissues, taken from the PaxDb protein-abundance database.',
    'tau': 'Tissue-specificity index (tau) adapted for protein abundance, computed for proteins present in at least five tissues; values near 0 indicate broad/ubiquitous expression and near 1 indicate tissue-specific expression.',
    'degrees': 'Node degree per residue in the Residue Interaction Network (RIN): the number of residue contacts (edges) of each residue.',
    'eccentricity': 'Eccentricity per residue in the RIN: the greatest shortest-path distance from the residue to any other residue.',
    'triangles': 'Number of triangles each residue participates in within the RIN (pairs of its contacting neighbours that also contact each other).',
    'eigenvalues': "Per-residue values from the spectral (eigenvalue) decomposition of the Residue Interaction Network's adjacency matrix.",
}
# e.g.  FEATURE_BRIEFS["Aggregation_8"] = "AGGRESCAN-derived aggregation propensity per residue."

ANALYSIS_SLUGS = {
    "btn-pca": "pca", "btn-umap": "umap", "btn-spider": "spider",
    "btn-heatmap": "heatmap", "btn-pin": "psn", "btn-violin": "violin",
    "btn-pgc": "pgc", "btn-profile": "profile", "btn-corr": "correlation",
}
ANALYSIS_BRIEFS = {
    "btn-pca": "Projects the standardized feature matrix onto uncorrelated principal components (the eigenvectors of the covariance matrix) that retain as much of the data's variance as possible, giving a 2D view of the selected proteins.",
    "btn-umap": "Uniform Manifold Approximation and Projection: a non-linear dimensionality-reduction method that captures high-dimensional structure in a 2D layout, yielding t-SNE-like results substantially faster.",
    "btn-spider": "Compares proteins against the human proteome or a user-defined set using per-feature Z-scores (deviation in standard deviations), highlighting the most divergent features while pruning within-category and highly-correlated (Spearman rho > 0.9) redundancy.",
    "btn-heatmap": "Standardized feature values per protein shown as a heatmap with hierarchical (Ward's method) clustering on Euclidean distance, grouping proteins and features with similar profiles.",
    "btn-pin": "A graph in which selected proteins are nodes and edges connect pairs with cosine similarity above 0.5, comparing proteins by the angle between their centered, preprocessed feature vectors rather than their magnitude.",
    "btn-violin": "Shows a feature's distribution via kernel density estimation, combining a box-plot-style summary (median, quartiles, IQR, outliers) with the estimated density curve.",
    "btn-pgc": "Ranks the features that most distinguish two protein groups by effect size (Cohen's d, with Hedges'/Durlak small-sample correction), using relative difference as a tiebreaker and reporting the top features.",
    "btn-profile": "Scatter plots of per-residue (min-max scaled) features along the sequence, letting the user inspect and compare specific regions of interest across proteins.",
    "btn-corr": "A pairwise Spearman's rho heatmap of the scalar features, giving an indicative view of monotonic correlations between features as a first-pass redundancy check (not a full multicollinearity analysis).",
}

def build_feature_descriptions(feature_ids):
    """feature_ids: ids actually present (numeric_df.columns + per-residue)."""
    ids = list(dict.fromkeys(list(feature_ids) + ALL_FEATURE_IDS))  # union, ordered
    return {fid: {"brief": FEATURE_BRIEFS.get(fid, ""), "url": feature_url(fid)} for fid in ids}

def build_analysis_descriptions():
    return {bid: {"brief": ANALYSIS_BRIEFS.get(bid, ""), "url": analysis_url(slug)}
            for bid, slug in ANALYSIS_SLUGS.items()}