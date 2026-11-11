# The PHANES database

*In Greek mythology, Phanes was the primordial god of creation. Etymologically, his name translates to “bring light”, as we hope to bring to your explorations of the human proteome.*

Here we introduce **PHANES (Proteome-wide Human Atlasof iNterpretable fEatureS)**, a comprehensive resource providing biologically interpretable protein descriptors for the human proteome. PHANES was designed to support exploratory analysis, comparative studies and machine-learning workflows by offering a unified and reusable representation of protein properties derived from sequence, structure, topology, stability and expression context.

For each human protein, PHANES integrates a collection of **sequence-derived physicochemical features**, including properties related to hydrophobicity, charge, aggregation propensity, nucleic acid binding, disorder, secondary-structure propensities and transmembrane regions. These descriptors capture intrinsic properties of amino-acid sequences and are provided both at residue level and as protein-level summaries.

To incorporate structural information, PHANES leverages predicted three-dimensional models from AlphaFoldDB. From these structures, we compute **geometric and surface-based features** describing protein size, compactness, accessibility and secondary-structure content, as well as confidence metrics derived from AlphaFold prediction scores. Structural information is further summarized using **graph-theoretical descriptors** by representing proteins as residue interaction networks, enabling the characterization of topology, centrality, hubness and clustering properties.

In addition to sequence and structure, PHANES includes **thermostability-related descriptors** derived from interaction energy statistics, providing an estimate of protein stability based on deviations from reference energy distributions. To place proteins in their biological context, PHANES also incorporates **tissue expression features**, capturing both overall abundance and tissue specificity across a wide range of human organs using large-scale proteomic datasets.

All features are precomputed for the human proteome and are accessible through an interactive web interface. Users can query proteins by UniProt accession, select and filter features of interest, visualize distributions and download structured feature matrices for downstream analysis. PHANES further allows users to upload protein structures in PDB format and compute the same set of descriptors on demand, enabling analysis of non-human, engineered or user-defined proteins (excluding expression-based features where appropriate).

PHANES is designed as a **general-purpose feature resource** rather than a predictor of specific biological phenomena. By prioritizing interpretability, modularity and completeness, it aims to complement modern embedding-based approaches and to facilitate transparent protein analysis, dataset construction and integration with statistical and machine-learning methods.

