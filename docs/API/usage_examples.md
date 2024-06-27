# Short example-drived guide to iLund4u API.  

iLund4u has a simple API allowing it programmatic usage from within a python program. Below we describe several Python snippets that mimic results of command-line calls.

**See detailed description of each class and method in the "Library" section.**


## Hotspot annotation mode



```python
import ilund4u

# Creating a parameters object and loading config
parameters = ilund4u.manager.Parameters()
parameters.load_config("standard")

# Example of changing a particular parameter
parameters.args["output_dir"] = "API_output_example"

# To turn off progress messages
parameters.args["verbose"] = False

# Creating a proteomes object and loading gff files  
proteomes = ilund4u.data_processing.Proteomes(parameters=parameters)

# Loading folder with gff files (for ex. the example-driven guide)
gff_folder = "ilund4u_data/guide/gff_files" 
proteomes.load_sequences_from_extended_gff(input_f=gff_folder)

# Running mmseqs on all encoded proteins and processing results
mmseqs_clusters = proteomes.mmseqs_cluster()
cluster_to_sequence = proteomes.process_mmseqs_results(mmseqs_clusters)

# Building proteome network, finding communities and defining protein classes        
proteome_network = proteomes.build_proteome_network(cluster_to_sequence)
proteomes.find_proteome_communities(proteome_network)
proteomes.define_protein_classes()

# Annotation of variable islands as regions with at least one non-conserved protein
proteomes.annotate_variable_islands()

# Building island network and finding hotspots there 
island_networks = proteomes.build_islands_network()
hotspots = proteomes.find_hotspots(island_networks)

# Additional functional annotation of proetins
hotspots.pyhmmer_annotation(proteomes)

# Merging hotspots by building hotspot network
hotspots.build_hotspot_network()

# Calculate statistics
hotspots.calculate_hotspot_statistics_and_get_annotation(proteomes)
hotspots.get_each_protein_group_statistics(proteomes)

# Saving results as database
database_manager = ilund4u.data_manager.DatabaseManager(parameters)
database_manager.build_database(proteomes, hotspots, "database_dirname")

# Visualisation with LoVis4u
drawing_manager = ilund4u.drawing.DrawingManager(proteomes, hotspots, parameters)
drawing_manager.plot_hotspot_communities()
drawing_manager.plot_proteome_communities()
```

---

## Protein and proteome annotation modes


``` py
import ilund4u


# Creating a parameters object and loading config
parameters = ilund4u.manager.Parameters()
parameters.load_config("standard")

# Example of changing a particular parameter
parameters.args["output_dir"] = "API_output_example"

# Database loading
database_manager = ilund4u.data_manager.DatabaseManager(parameters)
database = database_manager.load_database(db_path="path_to_the_db")

# Protein annotation mode
fasta_file_path = "ilund4u_data/guide/RloC.fa" 
database.protein_search_mode(query_fasta=fasta_file_path, query_label="RloC")

# Proteome annotation mode
gff_file_path = "ilund4u_data/guide/NC_001895.1.gff" 
database.proteome_annotation_mode(query_gff=gff_file_path)
```  
