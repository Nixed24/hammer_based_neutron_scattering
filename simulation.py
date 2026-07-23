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
def matrix_vector_multiply(matrix, vector):
    vector = np.squeeze(vector)
    output_vector = np.zeros(3)
    for i, row in enumerate(matrix):
        output_vector[i] = np.sum(row * vector)
    return output_vector
def do_rotation(vector, axis, _angle):
    if _angle > np.pi:
        angle = 2 * np.pi - _angle
    else:
        angle = _angle
    #print(f"angle: {angle}")
    if axis == "x":
        rotation_matrix = np.array([[1,0,0], \
                                    [0, np.cos(angle), -1 * np.sin(angle)],\
                                    [0, np.sin(angle), np.cos(angle)]])
    elif axis == "y":
        rotation_matrix = np.array([[np.cos(angle), 0, np.sin(angle)],\
                                    [0, 1, 0],\
                                    [-1 * np.sin(angle), 0, np.cos(angle)]])
    elif axis == "z":
        rotation_matrix = np.array([[np.cos(angle), -1 * np.sin(angle), 0],\
                                    [np.sin(angle), np.cos(angle), 0],\
                                    [0,0,1]])
    else:
        print(f"do_rotation(): No such axis '{axis}'")
        rotation_matrix = np.identity(3)
    rotated_vector = matrix_vector_multiply(rotation_matrix, vector)
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
    count_array = []
    for row in array_2d:
        count_array.append(np.count_nonzero(row))
    return count_array
def resolve_solid(solid):
    #TODO: prism behaves strangely off-diagonal
    solid_type = "none"
    solid_characteristics = {}
    plane_points_array = np.zeros((64,3,3))
    normal_vector_array = np.zeros((64,1,3))
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
    plane_points_array = np.round(plane_points_array)
    #normal_vector_array = np.round(normal_vector_array)
    dot_product_mesh = np.zeros((1, side_count)) # it's more like a matrix
    dot_product_array = np.zeros(side_count)
    for side_num in range (side_count):
        for other_side_num in range (side_count):
            dot_product_array[other_side_num] = np.dot(np.squeeze(normal_vector_array[side_num]), np.squeeze(normal_vector_array[other_side_num]))
        dot_product_mesh = np.vstack((dot_product_mesh, dot_product_array))
    dot_product_mesh = dot_product_mesh[1:] 
    """
        for row in range(dot_product_mesh.shape[0]):
            for col in range(dot_product_mesh.shape[1]):
                pass
    """
    dot_product_mesh = np.round(dot_product_mesh, 3)
    nonzero_list = count_nonzeroes(dot_product_mesh)
    print(dot_product_mesh)
    
    primary_side_numbers = []
    secondary_side_numbers = []
    if side_count == 5:
        solid_type = "prism"
        characteristic_side_number = nonzero_list.index(3)
        
        for i in range (len(dot_product_mesh[0])):
            if dot_product_mesh[characteristic_side_number][i] != 0 and dot_product_mesh[characteristic_side_number][i] != 1:
                # can also use sqrt(x) != x
                primary_side_numbers.append(i)
            if dot_product_mesh[characteristic_side_number][i] == 0:
                secondary_side_numbers.append(i)
                
        # Now do a basis transformation so that the prism fits nicely within XYZ (seperate rotations in 3 diff axes)
        # Condition: One primary side normal vector is a unit vector in x/y/z dir
        
        # all wrong
        # |
        # V
        """
        cos_in_x = np.sum(normal_vector_array[primary_side_numbers[0]] * [1,0,0])
        cos_in_y = np.sum(normal_vector_array[primary_side_numbers[0]] * [0,1,0]) #(it's just 90 - [x angle])
        cos_in_z = np.sum(normal_vector_array[primary_side_numbers[0]] * [0,0,1])
        theta_x = np.arccos(cos_in_x)
        theta_y = np.arccos(cos_in_y)
        theta_z = np.arccos(cos_in_z)
        """
        print(normal_vector_array[3][0])
        # It seems like it's MISSING a rotation?
        cos_xz = normal_vector_array[primary_side_numbers[0]][0][0] / np.sqrt(normal_vector_array[primary_side_numbers[0]][0][0]**2 + normal_vector_array[primary_side_numbers[0]][0][2]**2)
        cos_xy = normal_vector_array[primary_side_numbers[0]][0][0] / np.sqrt(normal_vector_array[primary_side_numbers[0]][0][0]**2 + normal_vector_array[primary_side_numbers[0]][0][1]**2)
        cos_yz = normal_vector_array[primary_side_numbers[0]][0][1] / np.sqrt(normal_vector_array[primary_side_numbers[0]][0][2]**2 + normal_vector_array[primary_side_numbers[0]][0][1]**2)
        theta_y = 0 # You only need two rotations to diagonalise, since after those two you'll be "into the page" so to speak on the XZ plane
        theta_z = np.arccos(cos_xy)
        theta_x = np.arccos(cos_yz)
        theta_giggle = 0 # in the current setup 0.9775 diagonalises it perfectly, where does this come from?
        print(f"{theta_x} | {theta_x * (180 / np.pi)}")
        print(f"{theta_y} | {theta_y * (180 / np.pi)}")
        print(f"{theta_z} | {theta_z * (180 / np.pi)}")
        
        diagonalised_normal_vector_array = np.copy(normal_vector_array)
        for i_vector, vector in enumerate(normal_vector_array):
            diagonalised_vector = do_rotation(vector, "x", theta_x)
            diagonalised_vector = do_rotation(diagonalised_vector, "y", theta_y)
            diagonalised_vector = do_rotation(diagonalised_vector, "z", theta_z)
            diagonalised_vector = do_rotation(diagonalised_vector, "z", theta_giggle)
            diagonalised_normal_vector_array[i_vector] = diagonalised_vector
        diagonalised_plane_points_array = np.zeros(plane_points_array.shape)
        for i_side, side in enumerate(plane_points_array):
            for i_point, point in enumerate(side):
                diagonalised_plane_point = do_rotation(point, "x", theta_x)
                diagonalised_plane_point = do_rotation(diagonalised_plane_point, "y", theta_y)
                #print(diagonalised_plane_point_xy)
                diagonalised_plane_point = do_rotation(diagonalised_plane_point, "z", theta_z)
                #print("B")
                #print(diagonalised_plane_point_xyz)
                diagonalised_plane_points_array[i_side][i_point] = diagonalised_plane_point
        print(primary_side_numbers)
        print(secondary_side_numbers)
        print(characteristic_side_number)
        prism_angle = np.arccos(np.sum(diagonalised_normal_vector_array[characteristic_side_number] * diagonalised_normal_vector_array[primary_side_numbers[0]]))
        if prism_angle > (np.pi * 0.5):
            prism_angle = np.pi - prism_angle
        diagonalised_plane_points_array = np.round(diagonalised_plane_points_array, 3)
        diagonalised_normal_vector_array = np.round(diagonalised_normal_vector_array, 3)
        # height = max Z - min Z
        # width = max Y - min Y
        # length = max X - min X
        # angle = acute arccos(diagonalised_normal_vectors[characteristic_side_number]) dotted with y unit vector
    if side_count == 6:
        pass
    if side_count > 6:
        pass
    return (prism_angle,\
            diagonalised_plane_points_array, plane_points_array,\
            diagonalised_normal_vector_array, normal_vector_array, dot_product_mesh)
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
        self.type = "cube"
        self.solid_parameter_dict = {"dimension", 1}
    def get_type(self):
        return
    def get_solid_parameters(self):
        return
    def set_type(self):
        #keep radius
        #discard most other things
        return
    def set_solid_parameter(self, parameter_name, parameter_value):
        return
    def set_solid_parameters(self, parameter_dict):
        return
class Solid:
    pass

test = resolve_solid(vmf_dict["world&0"]["solid&4"])