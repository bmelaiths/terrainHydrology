from __future__ import division
import numpy as np
import math
import scipy.spatial.distance
from mpl_toolkits.mplot3d import Axes3D

def random_point_line(num_points = 1):
    x = np.random.random(num_points)
    return np.reshape(x, (num_points,1))

def random_point_square(num_points = 1):
    x = np.random.random(num_points)
    y = np.random.random(num_points)
    return np.dstack((x,y))[0]


# if we only compare it doesn't matter if it's squared
def min_dist_squared(points, point):
    diff = points - np.array([point])
    return np.min(np.einsum('ij,ij->i',diff,diff))

class PoissonGenerator:
    def __init__(self, repeatPattern, first_point_zero):
        self.first_point_zero = first_point_zero
        self.repeatPattern = repeatPattern
        self.num_perms = (3 ** 2) if self.repeatPattern else 1


        self.zero_point = [0,0]
        
        
        self.random_point = random_point_square


    def first_point(self):
        if self.first_point_zero == True:
            return np.array(self.zero_point)
        return self.random_point(1)[0]

    def find_next_point(self, current_points, iterations_per_point):
        best_dist = 0
        best_point = []
        random_points = self.random_point(iterations_per_point)
        for new_point in random_points:
            dist = min_dist_squared(current_points, new_point)
            if dist > best_dist:
                best_dist = dist
                best_point = new_point
        return best_point

    def permute_point(self, point):
        out_array = np.array(point,ndmin = 2)
        if not self.repeatPattern:
            return out_array

        else :            
            for y in range(-1,2):
                for x in range(-1,2):
                    if y != 0 or x != 0:
                        perm_point = point+[x,y]
                        out_array = np.append(out_array, np.array(perm_point,ndmin = 2), axis = 0 )


        return out_array

    def find_point_set(self, num_points, num_iter, iterations_per_point, rotations, progress_notification = None):
        best_point_set = []
        best_dist_avg = 0
        self.rotations = 1

        for i in range(num_iter):
            if progress_notification != None:
                progress_notification(i / num_iter)
            points = self.permute_point(self.first_point())

            for i in range(num_points-1):
                next_point = self.find_next_point(points, iterations_per_point)
                points = np.append(points, self.permute_point(next_point), axis = 0)

            current_set_dist = 0

            if rotations > 1:
                points_permuted = np.copy(points)
                for rotation in range(1, rotations):
                    rot_angle = rotation * math.pi * 2.0 / rotations
                    s, c = math.sin(rot_angle), math.cos(rot_angle)
                    rot_matrix = np.matrix([[c, -s], [s, c]])
                    points_permuted = np.append(points_permuted, np.array(np.dot(points, rot_matrix)), axis = 0)
                current_set_dist = np.min(scipy.spatial.distance.pdist(points_permuted))
            else:
                current_set_dist = np.min(scipy.spatial.distance.pdist(points))

            if current_set_dist > best_dist_avg:
                best_dist_avg = current_set_dist
                best_point_set = points
        return best_point_set[::self.num_perms,:]


    