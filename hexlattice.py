import numpy as np
import csv
from hexlattice_func import *

# Example usage
s = 1  # size of the hexagon

############# Lattice size
nhex=100
#h = 12  # number of lines horizontally
#v = 12  # number of lines vertically

########## fraction settings
fraction_edge1=0.8
fraction_edge2=0.3
fraction_dbonds = 0.5  # Adjust the fraction here

# Define boundary type
boundary_type = 'x'  # options: 'open', 'x', 'y', 'xy'
    
try:
    h=nhex
    v=nhex
    hexgrid_path = 'hexgrid.csv'
    hexgrid = generate_hexgrid(s, h, v, save_path=hexgrid_path)
    
    # Define box size
    box_size = [s * h * np.sqrt(3), s * v * 1.5]

    # Paths for saving and loading
    hexgrid_path = 'hexgrid.csv'
    dist_matrix_path = 'dist_matrix.csv'
    edge_list1_path = 'edge_list1.csv'
    edge_list2_path = 'edge_list2.csv'

    # Calculate the distance matrix and save it
    if os.path.exists(dist_matrix_path):
        dist_matrix = load_csv(dist_matrix_path)
    else:
        dist_matrix = euclidean_distance_matrix(hexgrid, boundary_type, box_size, save_path=dist_matrix_path)
    
    # Create the edge lists and save them
    if os.path.exists(edge_list1_path):
        edge_list1 = load_edge_list(edge_list1_path)
    else:
        edge_list1 = create_edge_list(dist_matrix, 1, 1e-6, fraction_edge1, save_path=edge_list1_path)

    if os.path.exists(edge_list2_path):
        edge_list2 = load_edge_list(edge_list2_path)
    else:
        edge_list2 = create_edge_list(dist_matrix, 2, 1e-6, fraction_edge2, save_path=edge_list2_path)

    # Plot the hexgrid with edges
    #plot_hexgrid(hexgrid, edge_list1, edge_list2, show_vertex_numbers=True)

    output_file = 'particle_positions_'+str(h)+'by'+str(v)+'_'+boundary_type+'.txt'
    adjacency_list_path = 'Adjacency_list_'+str(h)+'by'+str(v)+'_'+boundary_type+'.txt'

    hexgrid = load_csv(hexgrid_path)
    edge_list1 = load_edge_list(edge_list1_path)
    edge_list2 = load_edge_list(edge_list2_path)

    # Merge edge lists
    edge_list = edge_list1 + edge_list2

    adjacency_list = generate_adjacency_list(hexgrid, edge_list, fraction_dbonds)

    create_adjacency_list_file(adjacency_list, adjacency_list_path)
    create_output_files(hexgrid, edge_list, output_file)

except ValueError as e:
    print(e)


