import numpy as np
import matplotlib.pyplot as plt
import csv
import os

################### all functions required to run hexlattice.py #########################
############################ designed for sequential run ################################
################ vertex numbering, corresponding edge numbering start from  1 ###########

def save_csv(data, file_path):
    with open(file_path, 'w', newline='') as file:
        writer = csv.writer(file, delimiter=',')
        writer.writerows(data)

def load_csv(file_path):
    data = []
    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            data.append([float(i) for i in row])
    return np.array(data)

def load_edge_list(file_path):
    with open(file_path, 'r') as file:
        reader = csv.reader(file, delimiter=',')
        return [tuple(map(int, row)) for row in reader]

def generate_hexgrid(s, h, v, save_path=None):
    # Ensure h and v are even for a perfect periodic lattice
    if h % 2 != 0 or v % 2 != 0:
        raise ValueError("h and v should be even for a perfect periodic lattice.")
    
    # Check if the hexgrid file already exists
    if save_path and os.path.exists(save_path):
        hexgrid = load_csv(save_path)
        return hexgrid

    # Generate the coordinates for the hexagonal grid
    x_coords = []
    y_coords = []

    for row in range(v):
        for col in range(h):
            if (row % 2 == 0 and col % 3 == 1) or (row % 2 == 1 and col % 3 == 2):
                continue
            x = s * (col * np.sqrt(3) + (row % 2) * (np.sqrt(3) / 2))
            y = s * (row * 1.5)
            x_coords.append(x)
            y_coords.append(y)

    hexgrid = np.column_stack((x_coords, y_coords))

    # Save the hexgrid to a file
    if save_path:
        save_csv(hexgrid, save_path)

    return hexgrid

def plot_hexgrid(hexgrid, edge_list1, edge_list2, show_vertex_numbers=False):
    fig, ax = plt.subplots(figsize=(8, 8))

    # Scatter plot of the vertices
    x_values, y_values = hexgrid[:, 0], hexgrid[:, 1]
    ax.scatter(x_values, y_values)

    # Plot the first set of edges
    for edge in edge_list1:
        x1, y1 = hexgrid[edge[0]-1]
        x2, y2 = hexgrid[edge[1]-1]
        ax.plot([x1, x2], [y1, y2], color='blue')

    # Plot the second set of edges
    for edge in edge_list2:
        x1, y1 = hexgrid[edge[0]-1]
        x2, y2 = hexgrid[edge[1]-1]
        ax.plot([x1, x2], [y1, y2], color='red')

    if show_vertex_numbers:
        for i, (x, y) in enumerate(zip(x_values, y_values), start=1):
            ax.text(x, y, str(i), fontsize=8, ha='center', va='center')

    ax.set_aspect('equal', 'box')
    ax.axis('off')
    plt.show()

def euclidean_distance_matrix(points, boundary_type, box_size, save_path=None):
    num_points = len(points)
    dist_matrix = np.zeros((num_points, num_points))
    
    for i in range(num_points):
        for j in range(i + 1, num_points):
            dist = euclidean_distance(points[i], points[j], boundary_type, box_size)
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist

    if save_path:
        save_csv(dist_matrix, save_path)
    
    return dist_matrix

def euclidean_distance(p1, p2, boundary_type, box_size):
    delta = np.abs(p1 - p2)
    if boundary_type == 'x' or boundary_type == 'xy':
        delta[0] = min(delta[0], box_size[0] - delta[0])
    if boundary_type == 'y' or boundary_type == 'xy':
        delta[1] = min(delta[1], box_size[1] - delta[1])
    return np.sqrt((delta ** 2).sum())

def create_edge_list(dist_matrix, order=1, tolerance=1e-6, fraction_edge=1, save_path=None):
    num_points = len(dist_matrix)
    edge_list = []

    res = dist_matrix.flatten()
    sorted_distances = np.unique(res.round(decimals=6))
    
    for i in range(num_points):
        for j in range(i + 1, num_points):
            if order < len(sorted_distances) and np.abs(dist_matrix[i, j] - sorted_distances[order]) < tolerance:
                edge_list.append((i+1, j+1))
    
    # Skip some edges based on the given fraction
    if fraction_edge < 1.0:
        np.random.shuffle(edge_list)
        edge_list = edge_list[:int(len(edge_list) * fraction_edge)]
    
    if save_path:
        save_csv(edge_list, save_path)
    
    return edge_list

def load_edge_list(file_path):
    edge_list = []
    with open(file_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            edge_list.append((int(row[0]), int(row[1])))
    return edge_list

def generate_adjacency_matrix(edge_list, num_vertices):
    adj_matrix = np.zeros((num_vertices, num_vertices))
    for edge in edge_list:
        adj_matrix[edge[0]-1, edge[1]-1] = 1
        adj_matrix[edge[1]-1, edge[0]-1] = 1  # Undirected graph
    return adj_matrix

def create_output_files(hexgrid, edge_list, output_file):
    # Create adjacency matrix
    num_vertices = len(hexgrid)
    adj_matrix = generate_adjacency_matrix(edge_list, num_vertices)

    # Find coordination numbers
    coordination_numbers = {}
    for i in range(num_vertices):
        coordination_numbers[i] = int(np.sum(adj_matrix[i]))  # Counting ones in the row

    # Generate output data
    #output_data = [['id', 'x', 'y', 'z']]
    output_data = []
    for i, (x, y) in enumerate(hexgrid, start=1):
        output_data.append([i, x, y, coordination_numbers[i - 1]])

    # Save output data to file
    save_csv(output_data, output_file)

def generate_adjacency_list(hexgrid, edge_list, fraction):
    adjacency_list = []
    for edge in edge_list:
        v1, v2 = edge
        x1, y1 = hexgrid[v1-1]
        x2, y2 = hexgrid[v2-1]
        nx = x2 - x1
        ny = y2 - y1
        rij = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        normal_x = (x2 - x1) / rij
        normal_y = (y2 - y1) / rij
        m = 1 if np.random.rand() < fraction else 0  # Assign 1 with the given fraction
        adjacency_list.append([v1, v2, nx, ny, m])
    return adjacency_list

def create_adjacency_list_file(adjacency_list, file_path):
#    header = ['id1', 'id2', 'nx', 'ny', 'm']
#    adjacency_list.insert(0, header)
    save_csv(adjacency_list, file_path)