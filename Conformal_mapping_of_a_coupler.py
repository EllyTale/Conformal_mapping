import numpy as np
from scipy.constants import epsilon_0
epsilon=11.75


def points_coupler(elements):
    '''
    This function converts elements of a coupler into complex points
    '''
    Z_points = np.zeros(len(elements) + 1)
    Z_points[0] = 1
    for i in range(len(elements)):
        Z_points[i+1] = Z_points[i] + elements[i]
    return np.array(Z_points)


def create_numerator_and_denumerator_points(points):
    '''
    This function helps find numerator and denumerator points
    '''
    numerator_point_ids = [p for p in np.arange(2, len(points)-2, 2)]
    denumerator_point_ids = [p for p in range(len(points)) if p not in numerator_point_ids]

    numerator_points = np.asarray(points)[numerator_point_ids]
    denumerator_points = np.asarray(points)[denumerator_point_ids]
    return numerator_points, denumerator_points


def function_for_points(points):
    '''
    This function create lists of points
    '''
    result = []
    for i in range(0, len(points)-2, 2):
        result.append(np.roll(points, -i))
    return result


def create_limits_Q(points):
    '''
    This function create list of limits for Q matrix
    '''
    result = list([])
    for k in range(1, len(points)-2, 2):
        result.append( [points[k], points[k+1]])
    return result


def create_limits_Phi(points):
    '''
    This function create list of limits for Phi matrix
    '''
    result=list([])
    for k in range(0, len(points)-2, 2):
        result.append( [points[k], points[k+1]])
    return result

def check_numerator_and_denumerator(numerator_points, denumerator_points, limits):
    '''
    This function checks and compares limits with numerator and denumerator points and it reterns new lists of numerator and denumerator points
    '''
    for k in range(len(limits)):
        if limits[k] in numerator_points:
            numerator_points.extend([limits[k]])
        else:
            denumerator_points.remove(limits[k])
    return numerator_points, denumerator_points


def gauss_chebyshev(numerator_points, denumerator_points, limits, n=100):
    '''
    This function counts Gauss-Chebyshev integral
    '''
    x = np.cos((2*np.arange(n)+1)*np.pi/(2*n))*(limits[1]-limits[0])*0.5+np.mean(limits)
    #y = np.ones(x.shape, np.complex)
    y = np.ones(x.shape)
    #print('x', x)
    for p in numerator_points:
        #print('num', p)
        y *= np.sqrt(np.abs(x-p))
    for p in denumerator_points:
        #print('den', p)
        y /= np.sqrt(np.abs(x-p))
    #print('y', y)
    return np.sum(y)*np.pi/n





class ConformalMapping:
    def __init__(self, elements):
        self.elements = np.asarray(elements)

    def cl(self):
        points = points_coupler(self.elements)
        print('Initial points', points)

        shape_of_matrix = int((len(points) - 2)/2)

        Q_mat = np.zeros((shape_of_matrix, shape_of_matrix))
        Phi_mat = np.zeros((shape_of_matrix, shape_of_matrix))


        # This part creates Q and Phi martix

        for i in range(shape_of_matrix):
            counter = 0
            phi_reference = 0

            list_of_points = function_for_points(points)[i]
            numerator_points, denumerator_points = create_numerator_and_denumerator_points(list_of_points)

            list_numerator_points_Q = list(numerator_points)
            list_denumerator_points_Q = list(denumerator_points)

            list_numerator_points_Phi = list(numerator_points)
            list_denumerator_points_Phi = list(denumerator_points)

            # This part counts Q of ground

            list_numerator_points_of_ground = list_numerator_points_Q
            list_denumerator_points_of_ground = list_denumerator_points_Q
            limits_of_ground = [list_of_points[len(list_of_points) - 1], list_of_points[0]]
            if i!=0:
                list_numerator_points_of_ground, list_denumerator_points_of_ground = check_numerator_and_denumerator(list_numerator_points_of_ground, list_denumerator_points_of_ground, limits_of_ground)
                Q_ground = gauss_chebyshev(list_numerator_points_of_ground, list_denumerator_points_of_ground, limits_of_ground, n=100)


            # Create list of limits for martix

            Q_list_of_limits = create_limits_Q(list_of_points)
            Phi_list_of_limits = create_limits_Phi(list_of_points)

            for j in range(shape_of_matrix):
                limits_Q = Q_list_of_limits[j]
                limits_Phi = Phi_list_of_limits[j]

                list_numerator_points_Q = list(numerator_points)
                list_denumerator_points_Q = list(denumerator_points)

                list_numerator_points_Phi = list(numerator_points)
                list_denumerator_points_Phi = list(denumerator_points)


                list_numerator_points_Q, list_denumerator_points_Q = check_numerator_and_denumerator(list_numerator_points_Q, list_denumerator_points_Q, limits_Q)
                list_numerator_points_Phi, list_denumerator_points_Phi = check_numerator_and_denumerator(list_numerator_points_Phi, list_denumerator_points_Phi, limits_Phi)

                if limits_Q == [points[len(points) - 1], points[0]]:
                    id1 = j
                    id2 = i
                else:
                    Q_mat[j][i] = gauss_chebyshev(list_numerator_points_Q, list_denumerator_points_Q, limits_Q, n=100)
                    counter=counter+Q_mat[j][i]

                Phi_mat[j][i] = phi_reference + gauss_chebyshev(list_numerator_points_Phi, list_denumerator_points_Phi, limits_Phi, n=100)
                phi_reference = Phi_mat[j][i]


            if i!=0:
                Q_mat[id1][id2] = Q_ground - counter


        print('Q = ', Q_mat)
        print('Phi = ', Phi_mat)

        Phi_inv = np.linalg.inv(Phi_mat)
        print(print('Phi_inv = ', Phi_inv))

        C = Q_mat*Phi_inv*(epsilon + 1)*epsilon_0

        return C
