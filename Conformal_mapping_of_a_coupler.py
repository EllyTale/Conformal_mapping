import numpy as np
from scipy.constants import epsilon_0
epsilon=11.75

'''
This function converts elements of a coupler into complex points
'''
def points_coupler(elements):
    Z_points=np.zeros(len(elements)+1)
    Z_points[0]=1
    for i in range(len(elements)):
        Z_points[i+1]=Z_points[i]+elements[i]

    #return np.array(Z_points, np.complex)
    return np.array(Z_points)

'''
This function helps find numerator and denumerator points
'''
def create_numerator_and_denumerator_points(points):

    numerator_point_ids = [p for p in np.arange(2,len(points)-2,2)]
    denumerator_point_ids = [p for p in range(len(points)) if p not in numerator_point_ids]

    numerator_points = np.asarray(points)[numerator_point_ids]
    denumerator_points = np.asarray(points)[denumerator_point_ids]

    return numerator_points, denumerator_points

'''
This function create lists of points
'''
def function_for_points(points):
    result = []
    for i in range(0, len(points)-2, 2):
        result.append(np.roll(points, -i))
    return result

'''
This function create list of limits for Q matrix
'''
def create_limits_Q(points):
    result=list([])
    for k in range(1, len(points)-2, 2):
        result.append( [points[k], points[k+1]])
    return result

'''
This function create list of limits for Phi matrix
'''
def create_limits_Phi(points):
    result=list([])
    for k in range(0, len(points)-2, 2):
        result.append( [points[k], points[k+1]])
    return result



'''
This function checks and compares limits with numerator and denumerator points and it reterns new lists of numerator and denumerator points
'''
def check_numerator_and_denumerator(numerator_points, denumerator_points, limits):
    for k in range(len(limits)):
        if limits[k] in numerator_points:
            numerator_points.extend([limits[k]])
        else:
            denumerator_points.remove(limits[k])

    return numerator_points, denumerator_points




'''
This function counts Gauss-Chebyshev integral
'''
def gauss_chebyshev(numerator_points, denumerator_points, limits, n=100):
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
        self.elements=np.asarray(elements)

    def C_l(self):
        points=points_coupler(self.elements)
        print(points)
        #numerator_points, denumerator_points = create_numerator_and_denumerator_points(points)

        shape_of_matrix=int((len(points)-2)/2)

        Q_mat = np.zeros((shape_of_matrix, shape_of_matrix))
        Phi_mat = np.zeros((shape_of_matrix, shape_of_matrix))

        '''
        This part creates Q martix
        '''

        for i in range(shape_of_matrix):
            counter=0
            '''
            Create lists of points
            '''
            list_of_points=function_for_points(points)[i]
            numerator_points, denumerator_points = create_numerator_and_denumerator_points(list_of_points)

            list_numerator_points=list(numerator_points)
            list_denumerator_points=list(denumerator_points)

            '''
            This part counts Q of ground
            '''
            list_numerator_points_of_ground=list_numerator_points
            list_denumerator_points_of_ground=list_denumerator_points
            limits_of_ground=[list_of_points[len(list_of_points)-1], list_of_points[0]]
            if i!=0:
                list_numerator_points_of_ground, list_denumerator_points_of_ground=check_numerator_and_denumerator(list_numerator_points_of_ground, list_denumerator_points_of_ground, limits_of_ground)
                Q_ground=gauss_chebyshev(list_numerator_points_of_ground, list_denumerator_points_of_ground, limits_of_ground, n=100)
                print('Q ground', Q_ground)

            #print('numerator points before', list_numerator_points)
            #print('denumerator points before', list_denumerator_points)
            '''
            Create list of limits for Q martix
            '''
            Q_list_of_limits=create_limits_Q(list_of_points)

            #print('List of limits', Q_list_of_limits)

            for j in range(shape_of_matrix):
                limits=Q_list_of_limits[j]
                list_numerator_points=list(numerator_points)
                list_denumerator_points=list(denumerator_points)
                '''
                This part checks and compares limits with numerator and denumerator
                '''
                list_numerator_points, list_denumerator_points=check_numerator_and_denumerator(list_numerator_points, list_denumerator_points, limits)
                print('limits Q', limits)
                #print('numerator Q', list_numerator_points)
                #print('denumerator Q', list_denumerator_points)
                if limits == [points[len(points)-1], points[0]]:
                    id1 = j
                    id2 = i


                else:
                    Q_mat[j][i]=gauss_chebyshev(list_numerator_points, list_denumerator_points, limits, n=100)
                    counter=counter+Q_mat[j][i]
                    #print('Q_elements',Q_mat[j][i])

            if i!=0:
                Q_mat[id1][id2]=Q_ground-counter


        '''
        This part creates Phi martix
        '''

        for i in range(shape_of_matrix):
            phi_reference=0

            '''
            Create list of points
            '''
            list_of_points=function_for_points(points)[i]
            numerator_points, denumerator_points = create_numerator_and_denumerator_points(list_of_points)

            list_numerator_points=list(numerator_points)
            list_denumerator_points=list(denumerator_points)
            #print('numerator points before', list_numerator_points)
            #print('denumerator points before', list_denumerator_points)

            '''
            Create list of limits for Phi matrix
            '''
            Phi_list_of_limits=create_limits_Phi(list_of_points)

            #print('List of limits', Q_list_of_limits)


            for j in range(shape_of_matrix):
                limits=Phi_list_of_limits[j]
                list_numerator_points=list(numerator_points)
                list_denumerator_points=list(denumerator_points)
                '''
                This part checks and compares limits with numerator and denumerator
                '''
                list_numerator_points, list_denumerator_points=check_numerator_and_denumerator(list_numerator_points, list_denumerator_points, limits)
                print('limits Q', limits)
                #print('numerator Phi', list_numerator_points)
                #print('denumerator Phi ', list_denumerator_points)



                Phi_mat[j][i]=phi_reference+gauss_chebyshev(list_numerator_points, list_denumerator_points, limits, n=100)
                phi_reference=Phi_mat[j][i]




        print('Q', Q_mat)
        print('Phi', Phi_mat)
        Phi_inv=np.linalg.inv(Phi_mat)
        print(print('Phi_inv', Phi_inv))

        C=Q_mat*Phi_inv*(epsilon+1)*epsilon_0



        return C
