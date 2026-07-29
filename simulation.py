#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 19 14:40:08 2026

@author: sainttux
"""

import vmf_deserialiser
import numpy as np

VMF_FILENAME = "test.vmf"

vmf_dict = vmf_deserialiser.deserialise_vmf(VMF_FILENAME)
def array_is_uniform(array):
    if array[0] == array[1]:
        for i, v in enumerate(array[:-1]):
            if v == array[i+1]:
                pass
            else:
                break
    else:
        return False
    return True
def get_fourth_vertex(plane_points):
    """
    Given a set of three vertices that correspond to a side that is perpendicular to 1 spatial (x/y/z) axis,
    find the fourth vertex of that side.

    Parameters
    ----------
    plane_points : np.ndarray
        Vertically stacked ndarray of points corresponding to the three vertices given in the plane field of the side

    Raises
    ------
    ValueError
        When the 'perpendicular to 1 axis' condition is not met.

    Returns
    -------
    fourth_point : np.ndarray
        A numpy array corresponding (x/y/z wise) to the 4th vertex of the side.

    """
    #Assuming plane is a rectangle w.r.t x-y-z
    #And the plane points are 3 vertices of that rectangle (they usually are in Hammer)
    #We could generalise this to work in any axis but we'd have to map our plane to this situation,
    #by choosing a different basis
    plane_points = np.array(plane_points)
    fourth_point = np.zeros(3)
    x_s, y_s, z_s = plane_points[:,0], plane_points[:,1], plane_points[:,2]
    point_array = np.array([x_s, y_s, z_s])
    for c_i, c in enumerate(point_array):
        if (np.average(c) - c[0]) == 0:
            const_axis_value = c[0]
            ax_nt = np.delete(np.arange(3), c_i) #axes_nontrivial
            points_of_interest = point_array
            fourth_point[c_i] = const_axis_value
    if not('ax_nt' in locals() and 'points_of_interest' in locals()):
        raise ValueError('get_fourth_vertex could not find constant axis. Are all your brushes cuboidal?')
    #We know the plane that the points form will have normal vector in the direction of axis specified by const_axis_value   
    dupe_table = np.zeros(points_of_interest.shape)
    for y in range(points_of_interest.shape[0]):
        for x in range(points_of_interest.shape[1]):
            # the value on a given axis with no dupes (i.e that value does not appear anywhere else on that axis)
            # will be the value on that axis of our missing point, but the constant axis will of course be constant.
            dupe_table[y][x] = np.sum(np.count_nonzero(points_of_interest[y] == points_of_interest[y][x]))
    for a in ax_nt:
        fourth_point[a] = point_array[a][dupe_table[a] == 1][0]
    return fourth_point
def interpolate_fourth_vertex(plane_points):
    x_s, y_s, z_s = plane_points[:,0], plane_points[:,1], plane_points[:,2]
    point_array = np.array([x_s, y_s, z_s])
    x_u, x_u_counts = np.unique(x_s, return_counts=True)
    fourth_vertex_x = x_u[np.argwhere(np.logical_or(x_u_counts == 1,x_u_counts == 3))]
    y_u, y_u_counts = np.unique(y_s, return_counts=True)
    fourth_vertex_y = y_u[np.argwhere(np.logical_or(y_u_counts == 1,y_u_counts == 3))]
    z_u, z_u_counts = np.unique(z_s, return_counts=True)
    fourth_vertex_z = z_u[np.argwhere(np.logical_or(z_u_counts == 1,z_u_counts == 3))]
    return np.squeeze(np.array([fourth_vertex_x, fourth_vertex_y, fourth_vertex_z]))
def do_rotation(vector, axis, _angle):
    if _angle > np.pi:
        angle = 2 * np.pi - _angle
    else:
        angle = _angle
    #print(f"angle: {angle}")
    if axis == "x":
        rotation_matrix = np.array([[1,0,0], \
                                    [0, np.cos(angle), -1 * np.sin(angle)],\
                                    [0, np.sin(angle), np.cos(angle)]], dtype=np.float64)
    elif axis == "y":
        rotation_matrix = np.array([[np.cos(angle), 0, np.sin(angle)],\
                                    [0, 1, 0],\
                                    [-1 * np.sin(angle), 0, np.cos(angle)]], dtype=np.float64)
    elif axis == "z":
        rotation_matrix = np.array([[np.cos(angle), -1 * np.sin(angle), 0],\
                                    [np.sin(angle), np.cos(angle), 0],\
                                    [0,0,1]], dtype=np.float64)
    else:
        print(f"do_rotation(): No such axis '{axis}'")
        rotation_matrix = np.identity(3)
    rotated_vector = np.matmul(rotation_matrix, vector)
    return rotated_vector
def magnitude(vector):
    square_sum = 0
    for value in vector:
        square_sum += value**2
    magnitude = np.sqrt(square_sum)
    return magnitude
def get_plane_points(side):
    plane_raw = side["plane"]
    plane_points = np.array(plane_raw.replace(")", "").replace("(", "").split(" ")).reshape((3,3))
    plane_points = np.array(plane_points, dtype=np.float32)
    return plane_points
def get_normal_vector(plane_points):
    line1 = plane_points[0] - plane_points[1]
    line2 = plane_points[0] - plane_points[2]
    normal_vector = np.cross(line1, line2)
    normal_vector = normal_vector / magnitude(normal_vector)
    return normal_vector
def count_nonzeroes(array_2d):
    count_array = np.zeros(array_2d.shape[0])
    for i_row, row in enumerate(array_2d):
        count_array[i_row] = np.count_nonzero(row)
    return count_array
def resolve_solid(solid):
    solid_type = "none"
    solid_dict = {}
    plane_points_array = np.zeros((256,3,3))
    normal_vector_array = np.zeros((256,1,3))
    side_count = 0
    counter = 0
    for key, value in zip(list(solid.keys()), list(solid.values())):
        if key[:4] == "side":
            current_plane = solid[key]
            current_plane_points = get_plane_points(current_plane)
            current_normal_vector = get_normal_vector(current_plane_points)
            plane_points_array[counter] = current_plane_points
            normal_vector_array[counter] = current_normal_vector
            side_count += 1
            counter += 1
    plane_points_array = plane_points_array[:counter]
    normal_vector_array = normal_vector_array[:counter]
    normal_vector_array = np.squeeze(normal_vector_array)
    dot_product_mesh = np.zeros((1, side_count)) # it's more like a matrix
    dot_product_array = np.zeros(side_count)
    for side_num in range (side_count):
        for other_side_num in range (side_count):
            dot_product_array[other_side_num] = np.dot(np.squeeze(normal_vector_array[side_num]), np.squeeze(normal_vector_array[other_side_num]))
        dot_product_mesh = np.vstack((dot_product_mesh, dot_product_array))
    dot_product_mesh = dot_product_mesh[1:] 
    dot_product_mesh_rounded = np.round(dot_product_mesh, 3)
    nonzero_list = count_nonzeroes(dot_product_mesh_rounded)
    if side_count == 5:
        primary_side_numbers = []
        secondary_side_numbers = []
        solid_type = "prism"
        characteristic_side_number = nonzero_list.tolist().index(3)
        
        for i in range (len(dot_product_mesh_rounded[0])):
            if dot_product_mesh_rounded[characteristic_side_number][i] != 0 and dot_product_mesh_rounded[characteristic_side_number][i] != 1:
                primary_side_numbers.append(i)
            if dot_product_mesh_rounded[characteristic_side_number][i] == 0:
                secondary_side_numbers.append(i)
        basis_vector_1 = normal_vector_array[primary_side_numbers[0]]
        basis_vector_2 = normal_vector_array[primary_side_numbers[1]]
        basis_vector_3 = normal_vector_array[secondary_side_numbers[0]]
        basis_matrix = np.transpose(np.vstack((basis_vector_1, basis_vector_2, basis_vector_3)))
        inverse = np.linalg.inv(basis_matrix)
        diagonalised_normal_vector_array = np.zeros(normal_vector_array.shape)
        for i_v, v in enumerate(normal_vector_array):
            new_v = np.matmul(inverse, v)
            diagonalised_normal_vector_array[i_v] = new_v
        diagonalised_plane_points_array = np.zeros((5,3,3))
        extended_diagonalised_plane_points_array = np.zeros((5,4,3))
        for i_side, side in enumerate(plane_points_array):
            for i_point, point in enumerate(side):
                diagonalised_point = np.matmul(inverse, point)
                diagonalised_plane_points_array[i_side][i_point] = diagonalised_point
                extended_diagonalised_plane_points_array[i_side][i_point] = diagonalised_point
        extended_diagonalised_plane_points_array = np.round(extended_diagonalised_plane_points_array, 1)
                #Diagonalised plane points will probably not be on integer units unfortunately
        
        side_centres = np.zeros((5,3))
        origin = np.zeros(3)
        for i_side, side in enumerate(extended_diagonalised_plane_points_array):
            fourth_point_coord = interpolate_fourth_vertex(side[:3])
            extended_diagonalised_plane_points_array[i_side][3] = fourth_point_coord
            for i in range (3):
                side_centres[i_side][i] = np.average(extended_diagonalised_plane_points_array[i_side][:,i])
        side_centres = np.delete(side_centres, characteristic_side_number, 0)
        side_centres = np.delete(side_centres, secondary_side_numbers, 0)
        # this way we know the origin will be at the centre of the characteristic side
        
        for i in range (3):
            origin[i] = np.average(side_centres[:,i])
        origin = np.matmul(origin, basis_matrix)
        
        diagonalised_normal_vector_array_rounded = np.round(diagonalised_normal_vector_array, 3)
        #WARNING : DO NOT USE TO COMPUTE INFORMATION ABOUT CHARACTERISTIC SIDE
        
        diagonalised_characteristic_side_parallel_axis = np.argmin(np.abs(diagonalised_normal_vector_array_rounded[characteristic_side_number]))
        parallel_points = diagonalised_plane_points_array[:,:,diagonalised_characteristic_side_parallel_axis]
        parallel_length = np.round(np.max(parallel_points) - np.min(parallel_points), 1)
        
        secondary_side_normal_vector_absolute = np.abs(diagonalised_normal_vector_array_rounded[secondary_side_numbers[0]])
        adjacent_secondary_side_parallel_axis = np.argwhere(secondary_side_normal_vector_absolute == np.min(secondary_side_normal_vector_absolute))[0]
        adjacent_secondary_side_parallel_unit_vector = np.zeros(3)
        adjacent_secondary_side_parallel_unit_vector[adjacent_secondary_side_parallel_axis] = 1
        adjacent_secondary_side_perp_axis = np.argwhere(secondary_side_normal_vector_absolute == np.min(secondary_side_normal_vector_absolute))[1]
        adjacent_secondary_side_perp_unit_vector = np.zeros(3)
        adjacent_secondary_side_perp_unit_vector[adjacent_secondary_side_perp_axis] = 1
        # Making sure our chosen parallel vector to get length will lead to the length of the adjacent side we do cos with (rather than 90 deg - cos)
        if np.isclose(np.sum(diagonalised_normal_vector_array_rounded[primary_side_numbers[0]] * adjacent_secondary_side_parallel_unit_vector), 0, atol=1e-4):
            pass
        else:
            adjacent_secondary_side_parallel_axis = np.argwhere(secondary_side_normal_vector_absolute == np.min(secondary_side_normal_vector_absolute))[1]
            adjacent_secondary_side_parallel_unit_vector = np.zeros(3)
            adjacent_secondary_side_parallel_unit_vector[adjacent_secondary_side_parallel_axis] = 1
            adjacent_secondary_side_perp_axis = np.argwhere(secondary_side_normal_vector_absolute == np.min(secondary_side_normal_vector_absolute))[0]
            adjacent_secondary_side_perp_unit_vector = np.zeros(3)
            adjacent_secondary_side_perp_unit_vector[adjacent_secondary_side_perp_axis] = 1
        adjacent_parallel_points = diagonalised_plane_points_array[:,:,adjacent_secondary_side_parallel_axis]
        adjacent_parallel_length = np.round(np.max(adjacent_parallel_points) - np.min(adjacent_parallel_points), 1)
        cos_angle = np.sum(diagonalised_normal_vector_array[characteristic_side_number]\
                           * diagonalised_normal_vector_array[primary_side_numbers[0]])
        if cos_angle < 0:
            cos_angle = -1 * cos_angle
        print(adjacent_parallel_length)
        characteristic_length = adjacent_parallel_length / (cos_angle)
        angle = np.arccos(cos_angle)
        solid_dict["type"] = solid_type
        solid_dict["angle"] = angle
        solid_dict["origin"] = origin # undiagonalised origin
        solid_dict["parallel_length"] = parallel_length
        solid_dict["characteristic_length"] = characteristic_length #diagonalised length
        #solid_dict["normals"] = diagonalised_normal_vector_array # un-diagonalised normals: prism normal SHOULD be the centre of the characteristic side in most cases
        solid_dict["basis_matrix"] = basis_matrix
        #solid_dict["characteristic_side"] = characteristic_side_number
        #solid_dict["primary_sides"] = primary_side_numbers
        #solid_dict["secondary_sides"] = secondary_side_numbers
        solid_dict["diagonalised_adjacent_parallel_dir"] = adjacent_secondary_side_parallel_unit_vector
        solid_dict["diagonalised_adjacent_perpendicular_dir"] = adjacent_secondary_side_perp_unit_vector
    if side_count == 6:
        i_allocated = []
        solid_type = "cuboid"
        # doesnt account for antiparallel and parallel being basically the same in our case
        dot_product_mesh_rounded = np.abs(dot_product_mesh_rounded)
        parallel_pairs = np.zeros((6,2), dtype=int)
        for i_row, row in enumerate(dot_product_mesh_rounded):
            parallel_pairs[i_row] = [np.where(row == 1)[0][0], np.where(row == 1)[0][1]]
        basis_vectors = []
        for side_number in np.unique(parallel_pairs[:,0]):
            basis_vectors.append(normal_vector_array[side_number])
        basis_matrix = np.transpose(np.vstack((basis_vectors[0], basis_vectors[1], basis_vectors[2])))
        inverse = np.linalg.inv(basis_matrix)
        diagonalised_normal_vector_array = np.zeros((6,3))
        for i_v, v in enumerate(normal_vector_array):
            new_v = np.matmul(inverse, v)
            diagonalised_normal_vector_array[i_v] = new_v
            
            
        extended_diagonalised_plane_points_array = np.zeros((6,4,3))
        diagonalised_plane_points_array = np.zeros((6,3,3))
        for i_side, side in enumerate(plane_points_array):
            for i_point, point in enumerate(side):
                diagonalised_point = np.matmul(inverse, point)
                diagonalised_plane_points_array[i_side][i_point] = diagonalised_point
                extended_diagonalised_plane_points_array[i_side][i_point] = diagonalised_point
        extended_diagonalised_plane_points_array = np.round(extended_diagonalised_plane_points_array, 1)
        side_centres = np.zeros((6,3))
        origin = np.zeros(3)
        for i_side, side in enumerate(extended_diagonalised_plane_points_array):
            fourth_point_coord = interpolate_fourth_vertex(side[:3])
            extended_diagonalised_plane_points_array[i_side][3] = fourth_point_coord
            for i in range (3):
                side_centres[i_side][i] = np.average(extended_diagonalised_plane_points_array[i_side][:,i])
        
        for i in range (3):
            origin[i] = np.average(side_centres[:,i])
        diagonalised_normal_vector_array_rounded = np.round(diagonalised_normal_vector_array)
        
        x_points = diagonalised_plane_points_array[:,:,0]
        y_points = diagonalised_plane_points_array[:,:,1]
        z_points = diagonalised_plane_points_array[:,:,2]
        
        x_len = np.round(np.max(x_points) - np.min(x_points), 1)
        y_len = np.round(np.max(y_points) - np.min(y_points), 1)
        z_len = np.round(np.max(z_points) - np.min(z_points), 1)
        
        x_sides = np.where(diagonalised_normal_vector_array_rounded[:,0] == 1) # or -1
        y_sides = np.where(diagonalised_normal_vector_array_rounded[:,1] == 1)
        z_sides = np.where(diagonalised_normal_vector_array_rounded[:,2] == 1)
        
        normdim_array = [[normal_vector_array[x_sides[0]], x_len],
                                  [normal_vector_array[y_sides[0]], y_len],
                                  [normal_vector_array[z_sides[0]], z_len]]
        solid_dict["type"] = solid_type
        solid_dict["origin"] = origin
        solid_dict["normdim"] = normdim_array
        solid_dict["basis_matrix"] = basis_matrix
    if side_count > 6: # not yet supported
        basis_vectors = []
        cap_side_numbers = np.argwhere(nonzero_list == 2)
        if len(cap_side_numbers) == 0:
            solid_type = "sphere"
            solid_dict["type"] = solid_type
        else:
            solid_type = "cylinder"
            basis_vectors.append(normal_vector_array[cap_side_numbers[0]][0])
            other_side_numbers = np.delete(np.arange(len(normal_vector_array)), cap_side_numbers)
            basis_vectors.append(normal_vector_array[other_side_numbers[0]])
            basis_vectors.append(np.cross(basis_vectors[0], basis_vectors[1]))
            basis_matrix = np.transpose(np.vstack((basis_vectors[0], basis_vectors[1], basis_vectors[2])))
            inverse = np.linalg.inv(basis_matrix)
            diagonalised_normal_vector_array = np.zeros(normal_vector_array.shape)
            for i_v, v in enumerate(normal_vector_array):
                new_v = np.matmul(inverse, v)
                diagonalised_normal_vector_array[i_v] = new_v
            cap_normal = np.round(diagonalised_normal_vector_array[cap_side_numbers[0]])[0]
            cap_axis = np.where(cap_normal != 0)[0][0]
            other_axes = np.delete(np.arange(3), cap_axis)
            diagonalised_plane_points_array = np.zeros(plane_points_array.shape)
            for i_side, side in enumerate(plane_points_array):
                for i_point, point in enumerate(side):
                    diagonalised_point = np.matmul(inverse, point)
                    diagonalised_plane_points_array[i_side][i_point] = diagonalised_point
            print(diagonalised_plane_points_array[:,:,cap_axis])
            length = np.round(np.max(diagonalised_plane_points_array[:,:,cap_axis]) - np.min(diagonalised_plane_points_array[:,:,cap_axis]), 1)
            diameter_estimate_1 = np.round(np.max(diagonalised_plane_points_array[:,:,other_axes[0]]) - np.min(diagonalised_plane_points_array[:,:,other_axes[0]]), 1)
            diameter_estimate_2 = np.round(np.max(diagonalised_plane_points_array[:,:,other_axes[1]]) - np.min(diagonalised_plane_points_array[:,:,other_axes[1]]), 1)
            solid_dict["type"] = solid_type
            solid_dict["length"] = length
            solid_dict["radius"] = np.average((diameter_estimate_1, diameter_estimate_2)) / 2 # this sucks
            # radius and origin how?
        #sphere and cylinder
    return solid_dict
        #dims
class Material:
    def __init__(self):
        self.atomic_number = 0
        self.microscopic_xsec = 0 #barn
        self.number_density = 0
        #macro = number_density * microscopic_xsec
        return
    def get_mean_free_path(self):
        return
class Solid_Base:
    def __init__(self):
        self.solid_dict = {}
        self.material = "vacuum"
        self.anti = False
    def attribute_exists(self, attribute):
        if attribute in self.solid_dict.keys():
            return True
        return False
    def get_attribute(self, attribute):
        if self.attribute_exists(attribute):
            return self.solid_dict[attribute]
        return None
    def get_solid_dict(self):
        return self.solid_dict
    def set_type(self):
        #keep radius
        #discard most other things
        return
    def merge_parameters(self, parameter_dict, overwrite=True):
        for key, value in zip(parameter_dict.keys(), parameter_dict.values()):
            if overwrite == False:
                if self.attribute_exists(key):
                    continue
                else:
                    self.solid_dict[key] = value
            else:
                self.solid_dict[key] = value
        return
    def point_is_inside(self, point):
        attribs = self.solid_dict
        inverse_basis_matrix = np.inverse(attribs["basis_matrix"])
        angle = attribs["angle"]
        if attribs["type"] == "prism":
            in_counter = 0
            max_len = np.max(attribs["dimensions"])
            diagonalised_point = np.matmul(point, inverse_basis_matrix)
            diagonalised_origin = np.matmul(attribs["origin"], inverse_basis_matrix)
            if magnitude(diagonalised_point - diagonalised_origin) > 2.1 * max_len:
                return False
            height = np.sin(angle) * attribs["characteristic_length"]
            length = np.cos(angle) * attribs["characteristic_length"]
            length_unit_vector = attribs["diagonalised_adjacent_parallel_dir"]
            length_axis = np.argmax(length_unit_vector)
            height_unit_vector = attribs["diagonalised_adjacent_perp_dir"]
            height_axis = np.argmax(height_unit_vector)
            depth_unit_vector = np.array([1,1,1]) - (np.array(length_axis) + np.array(height_axis))
            depth_axis = np.argmax(depth_unit_vector)
            depth = attribs["parallel_length"]
            depth_bounds = np.array([1,-1]) * depth / 2
            # we know for certain origin will be at midpoint in the parallel axis, it's at the centre of the characteristic side after all
            if np.abs(diagonalised_origin[depth_axis] - diagonalised_point[depth_axis]) < depth / 2:
                in_counter += 1
            # now our problem is 2D
            length_diff = diagonalised_origin[length_axis] - diagonalised_point[length_axis]
            height_diff = diagonalised_origin[height_axis] - diagonalised_point[height_axis]
            if not(np.abs(length_diff) < length / 2 and np.height_diff < height/2):
                return False
            #so it'd be in a cuboid made by shoving 2 of these prisms together
            diff_angle = np.arctan(height_diff / length_diff)
            if diff_angle <= angle:
                return True
            #TODO: implement epsilons?
        elif attribs["type"] == "cube":
            #TODO: this
            pass
        
            
                
class Solid:
    pass

test = resolve_solid(vmf_dict["world&0"]["solid&4"])