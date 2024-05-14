import numpy as np
import matplotlib.pyplot as plt
import random
import csv

def HexTile(s):
    return s * np.array([[np.sqrt(3)/2, 1/2] , [0, 1], [-np.sqrt(3)/2, 1/2], [-np.sqrt(3)/2, -1/2], [0, -1], [np.sqrt(3)/2, -1/2]])

def TranslateObject(p, translation):
    return [[point[0] + translation[0], point[1] + translation[1]] for point in p]

# regular open boundary hexagonal tiles
def HexGrid(s, h, v):
    hexgrid = []
    polygons = []
    for i in range(h):
        for j in range(v):
            translation = [s * (i * np.sqrt(3) + (j%2) * np.sqrt(3) / 2), 3 * s * j / 2]
            hexagon = HexTile(s) + translation
            # Check for overlapping vertices
            new_indices = []
            for vertex in hexagon:
                index = next((idx for idx, vec in enumerate(hexgrid) if np.allclose(vec, vertex)), None)
                if index is None:
                    hexgrid.append(vertex.tolist())
                    new_indices.append(len(hexgrid) - 1)
                else:
                    new_indices.append(index)
            # Append indices of non-overlapping vertices to polygons
            polygons.append(new_indices)
#    hlist = np.array(hexgrid)
    return hexgrid, polygons

def generate_edges(polygons, keep_fraction=1.0):
    edges = set()
    for polygon in polygons:
        for i in range(len(polygon)):
            edge = tuple(sorted([polygon[i], polygon[(i + 1) % len(polygon)]]))
            edges.add(edge)
    
    edges = list(edges)  # Remove duplicates

    if keep_fraction < 1.0:
        num_to_keep = int(len(edges) * keep_fraction)
        edges_to_keep = random.sample(edges, num_to_keep)
        edges = [edge for edge in edges if edge in edges_to_keep]
    return edges


def generate_crossing_edges(polygons, density=0.5):
    crossing_edges = set()
    for polygon in polygons:
        for i in range(len(polygon)):
            if random.random() < density:
                edge = tuple(sorted([polygon[i], polygon[(i + 2) % len(polygon)]]))
                crossing_edges.add(edge)
    crossing_edges=list(set(crossing_edges)) # remove duplicates
    return list(crossing_edges)

def find_boundary_points(hexgrid):
    boundary_indices = set()

    # Find the minimum and maximum x and y coordinates
    min_x = min(hexgrid, key=lambda point: point[0])[0]
    max_x = max(hexgrid, key=lambda point: point[0])[0]
    min_y = min(hexgrid, key=lambda point: point[1])[1]
    max_y = max(hexgrid, key=lambda point: point[1])[1]

    # Define the distance threshold for considering points on the boundary
    distance_threshold = 1  # Adjust as needed
    
    # Iterate over the hexgrid points and identify the boundary points
    for i, (x, y) in enumerate(hexgrid):
        if abs(x - min_x) < distance_threshold or abs(x - max_x) < distance_threshold \
           or abs(y - min_y) < distance_threshold or abs(y - max_y) < distance_threshold:
            boundary_indices.add(i)

    return list(boundary_indices)

# getting n-th digit of number
def get_digit(number, n):
    return number // 10**n % 10

# nh : horizontal cell numbers / nv : vertical cell numbers
def find_boundary_points_with_pairs(hexgrid, s,nh, nv,sel):
    boundary_indices = set()
    pairs = []

    # Find the minimum and maximum x and y coordinates
    min_x = min(hexgrid, key=lambda point: point[0])[0]
    max_x = max(hexgrid, key=lambda point: point[0])[0]
    min_y = min(hexgrid, key=lambda point: point[1])[1]
    max_y = max(hexgrid, key=lambda point: point[1])[1]

    # Define the distance threshold for considering points on the boundary
    distance_threshold = s*1  # Adjust as needed
    
    # Iterate over the hexgrid points and identify the boundary points
    for i, (x, y) in enumerate(hexgrid):
        if abs(x - min_x) < distance_threshold or abs(x - max_x) < distance_threshold \
           or abs(y - min_y) < distance_threshold or abs(y - max_y) < distance_threshold:
            boundary_indices.add(i)
    
    # tolerance for same x,y positions
    tolerance=0.01

    # Find pairs based on specified distances
    for i in boundary_indices:
        x, y = hexgrid[i]
        for j in boundary_indices:
            if j != i:
                x2, y2 = hexgrid[j]
                if (get_digit(sel, 1)==1):
                    # Find pairs for x boundaries
                    if (abs(abs(x - x2) - s* nh * np.sqrt(3)) < tolerance and  abs(y - y2) < tolerance):
                        if (x<x2):
                            pairs.append((i, j))
                        else:
                            pairs.append((j, i)) 
                if (get_digit(sel, 0)==1):
                    # Find pairs for y boundaries (even nv)
                    if (nv%2==0 and abs(y - y2) == s* nv * 1.5 and abs(x - x2) <tolerance):                    
                        if (y<y2):
                            pairs.append((i, j))
                        else:
                            pairs.append((j, i)) 
                    # Find pairs for y boundaries (odd nv)
                    elif (nv%2==1 and abs(y - y2) == s* (nv-1) * 1.5+s* 1 and abs(x - x2) <tolerance):                    
                        if (y<y2):
                            pairs.append((i, j))
                        else:
                            pairs.append((j, i)) 
                    elif (nv%2==1 and abs(y - y2) > s* (nv-1) * 1.5+s* 1 and abs(y - y2) < s* (nv-1)*1.5+s* 2+tolerance and abs(x - x2) <0.01):                    
                        if (y<y2):
                            pairs.append((i, j))
                        else:
                            pairs.append((j, i)) 
    pairs=list(set(pairs))
    return list(boundary_indices), pairs

def update_vertices_and_polygons(vertex_points, polygons, pairs):
    # Create a dictionary to map old vertex indices to new ones
    lenv=len(vertex_points)
    index_mapping = [[] for x in range(lenv)] # dictionary

    # Create a dictionary to map old and new vertex points
    updated_vertex_points = []

    # Remove boundary points from vertex_points
    for i, point in enumerate(vertex_points):
        if i not in [pair[0] for pair in pairs]:
            updated_vertex_points.append(vertex_points[i])

    # Add boundary points from pairs to vertex_points and create mappings (left: old / right: new address)
    for i, point in enumerate(vertex_points):
        index_mapping[i]=[i,False] # initial values, for points will be deleted.
        for j, point in enumerate(updated_vertex_points):
            if (vertex_points[i] == updated_vertex_points[j]):
                index_mapping[i]=[i,j]
    print(index_mapping)

    # Update polygons
    updated_polygons = polygons
    nhex=6
    npair=len(pairs)
    np=len(polygons)
    for polygon in range(0,np):
        for index in range(0,nhex):
            pvnum=polygons[polygon][index]
    # update polygon index by dictionary
    # update polygon list by pair (replacing deleted boundary to pair boundary)
            if (index_mapping[pvnum][1]==False):
                for i in range(0,npair):
                    # replace left to right one on the boundary
                    if (pvnum == pairs[i][0]):
                        npvnum = pairs[i][1]
                        if (index_mapping[npvnum][1]==False): # for periodic in both xy
                            for i in range(0,npair):
                                if (npvnum == pairs[i][0]):
                                    npvnum2 = pairs[i][1]
                            updated_polygons[polygon][index]=index_mapping[npvnum2][1]
                        else:
                            updated_polygons[polygon][index]=index_mapping[npvnum][1]
            else:
                updated_polygons[polygon][index]=index_mapping[pvnum][1]

    return updated_vertex_points, updated_polygons

# if vertex number starting from 1
def adjust_indices(edges):
    return [(u + 1, v + 1,w,x,y) for u, v,w,x,y in edges]

###### adding third column and vertex number (coordination number/neighbors)
def edges_to_vertex_count(edges, vertex):
    dim1 = len(vertex)
    vertex_count = [0] * dim1  # Initialize vertex count list with zeros
    for edge in edges:
        for v in edge[:2]:  # Limiting to only the first two components of each edge
            if 1 <= v <= dim1:  # Check if vertex number is within range
                vertex_count[v - 1] += 1  # Adjusting index for vertex number starting from 1
            else:
                print(f"Warning: Vertex number {v} is out of range")

    # Append counts to vertex coordinates
    vertex_with_count = [[i+1] + v + [vertex_count[i]] for i, v in enumerate(vertex)]
    return vertex_with_count


def add_random_weight(new_edges, probability, vertex_positions):
    weighted_edges = [(new_edge[0], new_edge[1], random.choice([0, 1])) for new_edge in new_edges]
    count_ones = sum(new_edge[2] for new_edge in weighted_edges)
    target_count_ones = round(len(weighted_edges) * probability)

    while count_ones != target_count_ones:
        new_edge = random.choice(weighted_edges)
        if new_edge[2] == 1:
            weighted_edges.remove(new_edge)
            weighted_edges.append((new_edge[0], new_edge[1], 0))
            count_ones -= 1
        else:
            weighted_edges.remove(new_edge)
            weighted_edges.append((new_edge[0], new_edge[1], 1))
            count_ones += 1

    # Calculate normal vector components (x and y direction) for each edge pair
    for i, edge in enumerate(weighted_edges):
        v1_x, v1_y = vertex_positions[edge[0]]
        v2_x, v2_y = vertex_positions[edge[1]]
        rij=np.sqrt((v1_x-v2_x)**2+(v1_y-v2_y)**2)
        normal_x = (v2_x-v1_x)/rij  # Normal vector component in x direction
        normal_y = (v2_y-v1_y)/rij  # Normal vector component in y direction
        weighted_edges[i] = (edge[0], edge[1],  normal_x, normal_y,edge[2])

    return weighted_edges

    
def plot_hexgrid_with_indices(hexgrid, polygons, boundary_indices, pairs):
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Plot polygons
    for polygon_indices in polygons:
        vertices = [idx for idx in polygon_indices if idx in hexgrid]  # Retrieve valid vertices for the polygon
        x_values = [hexgrid[idx][0] for idx in vertices]
        y_values = [hexgrid[idx][1] for idx in vertices]
        ax.fill(x_values, y_values, closed=True, edgecolor='black')
    
    # Plot boundary points with vertex indices
    for idx in boundary_indices:
        x, y = hexgrid[idx]
        ax.scatter(x, y, color='red')
        ax.text(x, y, str(idx), fontsize=8, ha='center', va='center')
    
    # Plot pairs
    for pair in pairs:
        idx1, idx2 = pair
        x1, y1 = hexgrid[idx1]
        x2, y2 = hexgrid[idx2]
        ax.plot([x1, x2], [y1, y2], color='blue')
    
    ax.set_aspect('equal', 'box')
    ax.axis('off')
    plt.show()

def plot_hexgrid(hexgrid, polygons, crossing_edges, boundary_points=None):
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Plot polygons
    for polygon_indices in polygons:
        vertices = [hexgrid[idx] for idx in polygon_indices if idx < len(hexgrid)]
        ax.fill(*zip(*vertices), closed=True, edgecolor='none')
    
    # Plot crossing edges
    for edge in crossing_edges:
        x_values = [hexgrid[edge[0]][0], hexgrid[edge[1]][0]]
        y_values = [hexgrid[edge[0]][1], hexgrid[edge[1]][1]]
        ax.plot(x_values, y_values, color='black')
    
    # Plot boundary points if provided
    if boundary_points:
        for index in boundary_points:
            x, y = hexgrid[index]
            ax.plot(x, y, marker='o', markersize=5, color='red')
    
    ax.set_aspect('equal', 'box')
    ax.axis('off')
    plt.show()

def plot_hexgrid_combined(hexgrid, polygons, crossing_edges, boundary_indices=None, pairs=None,hexgrid2=None):
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Plot polygons
    for polygon_indices in polygons:
        vertices = [hexgrid[idx] for idx in polygon_indices if idx < len(hexgrid)]
        ax.fill(*zip(*vertices), closed=False, edgecolor='none')
    
    # Plot crossing edges
    for edge in crossing_edges:
        x_values = [hexgrid[edge[0]][0], hexgrid[edge[1]][0]]
        y_values = [hexgrid[edge[0]][1], hexgrid[edge[1]][1]]
        ax.plot(x_values, y_values, color='black')
    
    # Plot boundary indices if provided
    if boundary_indices and hexgrid2:
        for idx, coord in zip(boundary_indices, hexgrid2):
            if idx < len(hexgrid2):
                x, y = hexgrid2[idx]
                ax.scatter(x, y, color='red')
                ax.text(x, y, str(idx), fontsize=8, ha='center', va='center')
            else:
                print(f"Warning: Boundary index {idx} is out of range")
    
    # Generate pairs coordinates from hexgrid if pairs are provided
    if pairs:
        pairs_coords = [[[hexgrid2[pair[0]][0], hexgrid2[pair[0]][1]],
                         [hexgrid2[pair[1]][0], hexgrid2[pair[1]][1]]] for pair in pairs]
        for pair, coords in zip(pairs, pairs_coords):
            idx1, idx2 = pair
            if idx1 < len(hexgrid2) and idx2 < len(hexgrid2):
                ax.plot([coords[0][0], coords[1][0]], [coords[0][1], coords[1][1]], color='blue')
            else:
                print(f"Warning: Pair indices {pair} contain index out of range")
      
    ax.set_aspect('equal', 'box')
    ax.axis('off')
    plt.show()
    
# Example usage:
s = 1  # Size of the hexagon (length of edge)
x_height = 20  # number of hexagons in x dir
y_height = 20  # number of hexagons in y dir 

hexgrid, polygons = HexGrid(s, x_height, y_height)

# Example usage:
sel = 0
# sel has two digit binary 0 (00) - open bd / 10 periodic in x / 1 (01) periodic in y / 11 periodic in both dir.
boundary_points, pairs  = find_boundary_points_with_pairs(hexgrid,s, x_height, y_height, sel)
#print("Boundary Points:", boundary_points)
print("Boundary Pairs:", pairs)
print("Polygons:", polygons)

edges = generate_edges(polygons,keep_fraction=0.8)
crossing_edges= generate_crossing_edges(polygons, density=0.3)
new_edges=edges+crossing_edges 
plot_hexgrid(hexgrid, polygons,new_edges,boundary_points)

# Test the function
updated_vertex_points, updated_polygons = update_vertices_and_polygons(hexgrid, polygons, pairs)
#print("Updated Vertex Points:", updated_vertex_points)
print("Updated Polygons:", updated_polygons)
edges = generate_edges(updated_polygons,keep_fraction=0.8)
crossing_edges= generate_crossing_edges(updated_polygons,density=0.3)
new_edges=edges+crossing_edges 
print(new_edges)

# plot removed vertex and updated at the same time
plot_hexgrid_combined(updated_vertex_points, updated_polygons,new_edges, boundary_points, pairs,hexgrid)
adjusted_edges=add_random_weight(new_edges,0.5,updated_vertex_points)
adjusted_edges = adjust_indices(adjusted_edges)
vertex_with_count = edges_to_vertex_count(adjusted_edges,updated_vertex_points)

with open(r'particle_positions_20by20_open.txt', 'w', newline="\n") as f:
    wr = csv.writer(f)
    wr.writerows(vertex_with_count)

    
with open(r'Adjacency_list_20by20_open.txt', 'w', newline="\n") as f:
    wr = csv.writer(f)
    wr.writerows(adjusted_edges)
print("save")
