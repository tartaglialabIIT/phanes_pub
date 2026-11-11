import os
from catgranuleFunctions import *
import concurrent.futures as cf
import networkx as nx
import argparse
from glob import glob
import os
import joblib
import json
from time import time
from column_mapping import columns_dic


# This calculator is now one stage of a longer pipeline (thermostability runs
# after it), so its 0-100 is remapped into a slice of the overall bar. Defaults
# to the full range, which is the previous behaviour.
PROGRESS_START = 0.0
PROGRESS_END = 100.0


def report(output_dir, pct, status):
    '''Write the current progress to progress.json in output_dir.
    The Dash app polls this file to drive the progress bar.'''
    try:
        pct = PROGRESS_START + (PROGRESS_END - PROGRESS_START) * pct / 100.0
        with open(os.path.join(output_dir, "progress.json"), "w") as f:
            json.dump({"progress": int(pct), "status": status}, f)
    except Exception:
        # never let progress reporting break the calculation
        pass


def get_nearest_atomic_distance(residue1, residue2):
    '''Calculates the minimum atomic distance of two residues'''
    min_dist = 14
    for atom1 in residue1:
        for atom2 in residue2:
            dist = np.linalg.norm(atom1.get_coord() - atom2.get_coord())
            if dist <= 12:
                if dist < min_dist:
                    min_dist = dist
    return min_dist


def build_graph_eduardo(G, pdb_structure, edge_distance=12):
    '''Build a Graph as in the https://dx.doi.org/10.1021/acs.jcim.9b01027
    The Node is the residue.
    An edge exists when the nearest atomic distance of two residues is <= 12 angstrom'''
    # Add nodes (residues) to the graph
    for model in pdb_structure:
        for chain in model:
            for residue in chain:
                G.add_node(residue)
    # Add edges to the graph
    for n in G.nodes:
        for n2 in G.nodes:
            if n != n2:
                if get_nearest_atomic_distance(n, n2) <= edge_distance:
                    G.add_edge(n, n2)

    return G


def compute_graph_theory_indicators_average_numpy(pdb_file):
    # Parse the PDB file
    print(pdb_file)
    name = os.path.splitext(os.path.basename(pdb_file))[0]
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)

    # Create a graph
    G = nx.Graph()
    G = build_graph_eduardo(G, structure)

    betweenness_centrality_mean = np.mean(np.array(list(nx.betweenness_centrality(G).values())))
    eigenvector_centrality_mean = np.mean(np.array(list(nx.eigenvector_centrality_numpy(G).values())))
    closeness_centrality_mean = np.mean(np.array(list(nx.closeness_centrality(G).values())))
    katz_centrality_mean = np.mean(np.array(list(nx.katz_centrality_numpy(G).values())))
    harmonic_centrality_mean = np.mean(np.array(list(nx.harmonic_centrality(G).values())))
    pagerank_mean = np.mean(np.array(list(nx.pagerank(G).values())))
    Hub_Scores_mean = np.mean(np.array(list(nx.hits(G)[0].values())))
    Authority_Scores_mean = np.mean(np.array(list(nx.hits(G)[1].values())))
    Local_Clustering_Coefficients_mean = np.mean(np.array(list(nx.clustering(G).values())))

    # Local Clustering measures how connected a residue's neighbors are to each other. |Global Clustering provides an overall measure of how interconnected the network is
    local_clustering_coefficient_profile = np.array(list(nx.clustering(G).values()))
    # The degree of a node is defined as the number of edges connected to it, which indicates how many direct connections (or neighbors) a node has
    degrees = [d for n, d in G.degree()]
    # In simple terms, it reflects how many edges are incident to a node, indicating its connectivity
    degree_centrality_profile = np.array(list(nx.degree_centrality(G).values()))
    degree_centrality_mean = np.sum(np.array(list(nx.degree_centrality(G).values()))) / len(nx.degree_centrality(G))

    betweenness_centrality_profile = np.array(list(nx.betweenness_centrality(G).values()))
    eigenvector_centrality_profile = np.array(list(nx.eigenvector_centrality_numpy(G).values()))
    closeness_centrality_profile = np.array(list(nx.closeness_centrality(G).values()))
    katz_centrality_profile = np.array(list(nx.katz_centrality_numpy(G).values()))
    harmonic_centrality_profile = np.array(list(nx.harmonic_centrality(G).values()))
    pagerank_profile = np.array(list(nx.pagerank(G).values()))
    Hub_Scores_profile = np.array(list(nx.hits(G)[0].values()))
    Authority_Scores_profile = np.array(list(nx.hits(G)[1].values()))

    # The shortest path length between two nodes is the minimum number of edges that must be traversed to get from one node to the other
    avg_shortest_path_length = np.nan
    # the longest shortest path between any two nodes in the graph
    diameter = np.nan
    if nx.is_connected(G):
        avg_shortest_path_length = nx.average_shortest_path_length(G)
        diameter = nx.diameter(G)
    # the greatest distance from a given node to any other node in the graph
    eccentricity_profile = np.array(list(nx.eccentricity(G).values()))
    #  the number of triangles that each node is a part of in a graph. A triangle is defined as a set of three nodes where each node is connected to the other two by edges
    triangles_profile = np.array(list(nx.triangles(G).values()))
    # adjacency matrix of a graph provides significant insights into its structural properties, dynamics
    adj_matrix = nx.to_numpy_array(G)
    eigenvalues_profile = np.linalg.eigvals(adj_matrix)
    # entropy measures the level of uncertainty or disorder within the degree distribution of the nodes. It quantifies how uniformly the degrees of the nodes are distributed
    degrees_profile = [d for n, d in G.degree()]
    entropy = -sum(p * np.log2(p) for p in np.array(degrees) / np.sum(degrees) if p > 0)

    return [name, len(G.nodes()), len(G.edges()), nx.density(G), nx.number_connected_components(G),
            len(max(nx.connected_components(G), key=len)), sum(dict(G.degree()).values()) / len(G),
            nx.degree_assortativity_coefficient(G), betweenness_centrality_mean, eigenvector_centrality_mean,
            closeness_centrality_mean, katz_centrality_mean, harmonic_centrality_mean, pagerank_mean,
            Hub_Scores_mean, Authority_Scores_mean, Local_Clustering_Coefficients_mean,
            nx.average_clustering(G), nx.transitivity(G), local_clustering_coefficient_profile, degrees_profile,
            degree_centrality_profile, degree_centrality_mean,
            betweenness_centrality_profile, eigenvector_centrality_profile, closeness_centrality_profile,
            katz_centrality_profile,
            harmonic_centrality_profile, pagerank_profile, Hub_Scores_profile, Authority_Scores_profile,
            avg_shortest_path_length,
            diameter, eccentricity_profile, triangles_profile, eigenvalues_profile, entropy]


def build_graph_parallel(pdb_list, progress_cb=None):
    '''For a pdb list this calculates all the graph indicators as a data list
    This should either be a list of .pdb or pdb.gz files.
    If progress_cb is given, it is called as progress_cb(done, total) after each
    PDB finishes. Results may come back out of order, which is fine because the
    caller re-indexes by id before concatenating.'''
    data_par = []
    total = len(pdb_list)
    with cf.ProcessPoolExecutor(max_workers=6) as executor:
        futures = [executor.submit(compute_graph_theory_indicators_average_numpy, p)
                   for p in pdb_list]
        for done, fut in enumerate(cf.as_completed(futures), start=1):
            data_par.append(fut.result())
            if progress_cb is not None:
                progress_cb(done, total)
    return data_par


def build_graph_enriched_parallel(pdb_list):
    '''For a pdb list this calculates all the graph indicators as a data list
    This should either be a list of .pdb or pdb.gz files'''
    data_par = []
    with cf.ProcessPoolExecutor(max_workers=6) as executor:
        data_par = list(executor.map(compute_graph_theory_indicators_enriched, pdb_list))
    data_par.to_pickle('graph_enrichment_arr.pkl')
    return data_par


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="This program calculates physical/chemical and alphafold features")
    parser.add_argument('-D', "--pdb_dir", help="directory that contains pdb files", required=True)
    parser.add_argument('-O', "--output_dir", help="Path to the output directory", required=True)
    parser.add_argument("--progress-start", type=float, default=0.0,
                        help="map this script's 0-100 progress into [start, end]")
    parser.add_argument("--progress-end", type=float, default=100.0)

    # Parse the command-line arguments
    args = parser.parse_args()

    # module scope here, so this rebinds the globals report() reads
    PROGRESS_START, PROGRESS_END = args.progress_start, args.progress_end

    # Access the parsed arguments
    code_dir = os.getcwd()
    scales_dir = './src/ChemicalPhysicalScales_Py_dictionary'
    classifiers_dir = './src/TRAINED_MODELS/'
    af_code_dir = './src/AlphaFold/'
    pdb_dir = args.pdb_dir
    output_dir = args.output_dir


    pdb_files = glob(f"{pdb_dir.rstrip('/')}/*pdb")
    start = time()
    report(output_dir, 5, "Reading PDB files...")
    sequences, ids = get_sequences_from_pdbs(pdb_files)
    # get the dataframe with the physical/chemical properties
    report(output_dir, 20, "Computing physico-chemical properties...")
    # get_physical_chemical() instead of get_physical_chemical_properties():
    # both average the same per-residue matrix, but this one also returns that
    # matrix, giving us the per-residue profile of every scale alongside the
    # scalar mean. The scalar columns are unchanged -- compute_chemphys and
    # compute_chemphysProperties have identical bodies and differ only in what
    # they return. The last two arguments are unused by the function.
    pc_df = get_physical_chemical(sequences, ids, scales_dir, classifiers_dir, None)
    end_pc = time()
    print(f'pc = {end_pc - start} seconds')
    report(output_dir, 40, "Extracting AlphaFold features...")
    af_df = get_alphafold_features_from_file(af_code_dir, pdb_files)
    # add plddt columns
    af_complete_df = add_plddt_to_af_features_from_file(pdb_files, af_df)
    postprocessed_af_df = ProcessAF_Features_from_file(af_complete_df)
    # set the uniprot _ID as the index
    postprocessed_af_df.set_index('Uniprot_ID', inplace=True)
    end_af = time()
    print(f'af = {end_af - end_pc} seconds')
    graph_cols = ['id', 'number_of_nodes', 'number_of_edges', 'Density', 'number_of_connected_components',
                  'size_of_largest_connected_component', 'avgerage_degree',
                  'degree_assortativity', 'betweeness_centrality',
                  'eigenvector_centrality', 'closeness_centrality',
                  'katz_centrality', 'harmonic_centrality', 'pagerank', 'hub_scores',
                  'authority_scores', 'local_clustering_coefficient',
                  'average_clustering_coefficient', 'global_transitivity', 'local_clustering_coefficient_profile',
                  'degrees_profile', 'degree_centrality_profile', 'degree_centrality',
                  'betweeness_centrality_profile', 'eigenvector_centrality_profile', 'closeness_centrality_profile',
                  'katz_centrality_profile', 'harmonic_centrality_profile', 'pagerank_profile', 'hub_scores_profile',
                  'authority_scores_profile', 'average_shortest_path_length', 'diameter', 'eccentricity_profile',
                  'triangles_profile', 'eigenvalues_profile', 'entropy']
    # graph building is the slow phase -> report per-file progress over the 65-95 band
    GRAPH_START, GRAPH_END = 65, 95

    def graph_progress(done, total):
        pct = GRAPH_START + (GRAPH_END - GRAPH_START) * done / total
        report(output_dir, pct, f"Building protein graphs ({done}/{total})...")

    report(output_dir, GRAPH_START, "Building protein graphs...")
    data = build_graph_parallel(pdb_files, progress_cb=graph_progress)
    graphs_df = pd.DataFrame(data=data, columns=graph_cols)
    graphs_df.set_index('id', inplace=True)
    end_graph = time()
    print(f'af = {end_graph - end_af} seconds')
    report(output_dir, 97, "Combining results...")
    cmb_df = pd.concat([pc_df, postprocessed_af_df, graphs_df], axis=1)
    # rename raw scale-file column names to the short user-facing labels.
    # get_physical_chemical() suffixes its per-residue columns "_prof"; the app
    # and profet_numeric_v6_profiles.pkl both use "<label>_profile", which is
    # also what the graph profiles already use.
    def _label(col):
        if col.endswith("_prof"):
            stem = col[: -len("_prof")]
            return f"{columns_dic.get(stem, stem)}_profile"
        return columns_dic.get(col, col)

    cmb_df = cmb_df.rename(columns=_label)
    end = time()
    print(f'time elapsed is {(end - start) / 60} minutes')
    joblib.dump(cmb_df, f'{output_dir.rstrip("/")}/cmb_df.pkl')
    report(output_dir, 100, "Done")
