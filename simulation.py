        #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 19 14:40:08 2026

@author: sainttux
"""
from scipy import constants as constants # We use this library for Avogadro's constant
import vmf_deserialiser
import numpy as np

AVOGADRO_CONSTANT = constants.Avogadro
PLANCK_CONSTANT = constants.Planck # 6.63e-34 Js
NEUTRON_MASS = constants.neutron_mass # 1.6749e-27 kg
PROTON_MASS = constants.proton_mass # 1.6726e-27 kg
ELECTRON_MASS = constants.electron_mass # 9.11e-31 kg
ELECTRON_VOLT = constants.electron_volt # 1.6e-19 J

#VALID_PARTICLE_TYPES = ["neutron", "beta_minus", "beta_plus", "gamma"]
VALID_PARTICLE_TYPES = ["neutron"]

class Material_Profile:
    def __init__(self, name, atomic_number=None, atomic_mass=None, density=None, molar_mass=None):
        self.name = name
        self.density = density # density in grams per centimetre cubed
        self.Z = atomic_number # u
        self.A = atomic_mass # u
        self.molar_mass = molar_mass # molar mass in grams per mole
        self.atomic_radius = (1.2) * (self.A ** (1/3)) * (10 ** -15) #femtometres
        self.temperature = None
        self.number_density = ((self.density) / (self.molar_mass)) * (AVOGADRO_CONSTANT)
        return
    def read_in_material_properties(self, filename):
        return
    def get_mean_free_path(self, particle):
        if particle.type == "neutron":
            if self.temperature is None: #For liquids/gases only
                pass
            else:
                particle_momentum = particle.kinetic_energy / (2 * particle.mass)
                db_wavelength = PLANCK_CONSTANT / particle_momentum
                effective_radius = self.atomic_radius + (db_wavelength) / (2 * np.pi)
                microscopic_xsec = 2 * np.pi * (effective_radius**2)
                macroscopic_xsec = microscopic_xsec * self.number_density * 10**(-24) # cm^-1
                return (1 / macroscopic_xsec) * 0.01 # metres
        return

material_graphite = Material_Profile("graphite", 6, 12, 1.67, 12.011)
material_vacuum = Material_Profile("vacuum", 1e-10, 2e-10, 1e-30, 1e-60)
material_lead = Material_Profile("lead", 82, 207, 11.35, 207.2)
material_water = Material_Profile("water", 10, 34, 1, 18.0153)

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
def random_unit_sphere_points(size):
    phi_values = np.random.uniform(size = size) * 2 * np.pi
    theta_values = np.arccos(1-2*np.random.uniform(size = size)) 
    # np.random.uniform() won't be uniform in a spherical coords conversion, we need to put it through this function to uniformise it.
    x, y, z = np.sin(theta_values) * np.cos(phi_values), np.sin(theta_values) * np.sin(phi_values), np.cos(theta_values) #Converting to cartesian
    vectors = np.hstack((x,y,z))
    return vectors
def exp_distribution_float(mean, size):
    radii = -1 * mean * np.log(np.random.uniform(size=size))
    return radii
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
def resolve_solid(solid, origin=None):
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
        if key == "id":
            solid_dict["id"] = value
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
        if origin is None: # This sucks and I need to figure out how to do it properly, for now we'll just have
        #all brushes be brush ents and then get the origin from there.
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
            print(origin)
            origin = np.matmul(origin, basis_matrix) #undiagonalising origin
            print(origin)
        
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
        solid_dict["diagonalised_adjacent_perp_dir"] = adjacent_secondary_side_perp_unit_vector
    if side_count == 6:
        i_allocated = []
        solid_type = "cuboid"
        # doesnt account for antiparallel and parallel being basically the same in our case
        dot_product_mesh_rounded = np.abs(dot_product_mesh_rounded)
        parallel_pairs = np.zeros((6,2), dtype=int)
        print(dot_product_mesh_rounded)
        for i_row, row in enumerate(dot_product_mesh_rounded):
            parallel_pairs[i_row] = [np.argwhere(row == 1)[0][0], np.argwhere(row == 1)[1][0]]
        basis_vectors = []
        for side_number in np.unique(parallel_pairs[:,0]):
            basis_vectors.append(normal_vector_array[side_number])
        already_diagonal_check_array = np.array(np.abs(basis_vectors))
        if np.count_nonzero(already_diagonal_check_array) == 3: 
            # Only three non-zero elements in basis matrix -> just a 90 * N degree rotation on all axes -> no need to rotate
            basis_matrix = np.identity(3)
            inverse = basis_matrix
        else:
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
        if origin is None:
            side_centres = np.zeros((6,3))
            undiagonalised_side_centres = np.zeros((6,3))
            origin = np.zeros(3)
            for i_side, side in enumerate(extended_diagonalised_plane_points_array):
                fourth_point_coord = interpolate_fourth_vertex(side[:3])
                extended_diagonalised_plane_points_array[i_side][3] = fourth_point_coord
                for i in range (3):
                    side_centres[i_side][i] = np.average(extended_diagonalised_plane_points_array[i_side][:,i])
                undiagonalised_side_centres[i_side] = np.matmul(side_centres[i_side], basis_matrix)
            side_centres_rounded = np.round(side_centres, 1)
            print(side_centres_rounded)
            for i in range (3):
                origin[i] = np.average(undiagonalised_side_centres[:,i])
            origin = np.matmul(basis_matrix, origin)
            print(origin)
        
        diagonalised_normal_vector_array_rounded = np.round(diagonalised_normal_vector_array)
        
        
        solid_dict["type"] = solid_type
        solid_dict["origin"] = origin # undiagonalised
        solid_dict["d_side_centres"] = side_centres # diagonalised
        solid_dict["d_side_normals"] = diagonalised_normal_vector_array # diagonalised
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
class Particle:
    def __init__(self, name="neutron"):
        self.pos = np.array([0,0,0])
        self.rest_mass = 0
        self.energy = 0
        self.kinetic_energy = 0
        self.name = name
        if name in VALID_PARTICLE_TYPES:
            exec(f"init_{name}()")
        else:
            print(f"Particle not recognised: '{name}', defaulting to neutron")
            self.init_neutron()
    def init_neutron(self):
        self.name = "neutron"
        self.rest_mass = NEUTRON_MASS
        self.charge = 0
        return
    # get mfp formula function
    # get angular scattering cross-section (Klein-Nishina is an example of this)
    

class Solid_Base:
    def __init__(self):
        self.solid_dict = {}
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
        inverse_basis_matrix = np.linalg.inv(attribs["basis_matrix"])
        if attribs["type"] == "prism":
            angle = attribs["angle"]
            max_len = attribs["characteristic_length"]
            diagonalised_point = np.matmul(point, inverse_basis_matrix)
            diagonalised_origin = np.matmul(attribs["origin"], inverse_basis_matrix)
            if magnitude(diagonalised_point - diagonalised_origin) > 1.1 * max_len:
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
            # we know for certain origin will be at midpoint in the parallel axis, it's at the centre of the characteristic side after all
            if not (np.abs(diagonalised_origin[depth_axis] - diagonalised_point[depth_axis]) < depth / 2):
                return False
            # now our problem is 2D
            length_diff = diagonalised_origin[length_axis] - diagonalised_point[length_axis]
            height_diff = diagonalised_origin[height_axis] - diagonalised_point[height_axis]
            if not(np.abs(length_diff) < length / 2 and np.abs(height_diff) < height/2):
                return False
            #so it'd be in a cuboid made by shoving 2 of these prisms together
            diff_angle = np.arctan(height_diff / length_diff)
            if diff_angle <= angle:
                return True
            return False
        elif attribs["type"] == "cuboid":
            diagonalised_point = np.matmul(point, inverse_basis_matrix)
            diagonalised_origin = np.matmul(attribs["origin"], inverse_basis_matrix)
            side_centres = attribs["d_side_centres"]
            x_centre_points = attribs["d_side_centres"][:,0]
            y_centre_points = attribs["d_side_centres"][:,1]
            z_centre_points = attribs["d_side_centres"][:,2]
            x_len = np.max(x_centre_points) - np.min(x_centre_points)
            y_len = np.max(y_centre_points) - np.min(y_centre_points)
            z_len = np.max(z_centre_points) - np.min(z_centre_points)
            diff = np.abs(diagonalised_point - diagonalised_origin)
            if (diff[0] < (x_len/2) and diff[1] < (y_len/2) and diff[2] < (z_len/2)):
                return True
            return False
        #TODO: implement epsilons?
class Solid: #solid_base with physics implemented (Material basically)
    def __init__(self, solid, material):
        self.solid_profile = solid #Solid_Base
        self.material_profile = material #Material_Profile
        return
    def point_is_inside(self, point):
        return self.solid_profile.point_is_inside(point)
    def get_mfp(self, particle):
        return self.material_profile.get_mean_free_path(particle)
    pass
def solid_base_from_brush_ent(entity_dict):
    solid_dict = resolve_solid(entity_dict["solid&0"], np.array(entity_dict["origin"].split(" "), dtype=float))
    solid_instance = Solid_Base()
    solid_instance.merge_parameters(solid_dict)
    return solid_instance
def solid_from_world_brush():
    return
def simulation_construct_solids(): # Resolve each solid in the VMF into a Solid instance, using the texture of the first face to get the material
    return
def simulation_create_particle_stack(): 
    # Find info_targets (tentatively) with targetnames that match particle types -- these are our sources.
    # create the particles that originate from our sources at t_0 using the keyvalue indicating their intensity (name: amount)
    # Project these particles by their mfp from the source using the source's angles as the direction
    return
def simulation_main_loop():
    return
"""
test = resolve_solid(vmf_dict["world&0"]["solid&0"])
test_solid = Solid_Base()
test_solid.merge_parameters(test)
test_bool = test_solid.point_is_inside(np.array([-96,160,-32]))
"""
SIMULATION_BOUNDS = np.array([[-1000, 1000], [-1000, 1000], [-1000, 1000]])
test_solid_2 = solid_base_from_brush_ent(vmf_dict["entity&0"])
test_solid_2_dict = test_solid_2.get_solid_dict()
test_bool_2 = test_solid_2.point_is_inside(np.array([144,176,-176]))
