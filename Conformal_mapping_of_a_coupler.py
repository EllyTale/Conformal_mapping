import numpy as np
from scipy.constants import epsilon_0
epsilon=11.75

def create_numerator_and_denumerator_points(points):
    '''
    This function helps find numerator and denumerator points
    '''
    shape_of_matrix=(len(points)-2)/2
    #len(points)
    #print('Shape of matrix=', shape_of_matrix)
    create_numerator_points=np.zeros(int(shape_of_matrix)-1)
    create_denumerator_points=np.zeros(len(points)-(int(shape_of_matrix)-1))
    j=0
    if shape_of_matrix > 1:
        for i in range(len(create_numerator_points)):
            create_numerator_points[i]=points[j+2]
            j=j+2
    else:
        create_numerator_points=None

    if shape_of_matrix > 1:

        create_denumerator_points[0]=points[0]
        create_denumerator_points[1]=points[1]

        create_denumerator_points[len(create_denumerator_points)-1]=points[len(points)-1]
        create_denumerator_points[len(create_denumerator_points)-2]=points[len(points)-2]
        create_denumerator_points[len(create_denumerator_points)-3]=points[len(points)-3]

        k=3
        for i in range(2, len(create_denumerator_points)-3):
            create_denumerator_points[i]=points[k]
            k=k+2
    else:
         for i in range(len(create_denumerator_points)):
                create_denumerator_points[i]=points[i]

    numerator_points=create_numerator_points
    denumerator_points=create_denumerator_points

    #print('Numerator points', numerator_points)
    #print('Deumerator points', denumerator_points)

    return numerator_points, denumerator_points


def function_for_points(points):
    '''
    This function create lists of points
    '''
    shape_of_matrix=(len(points)-2)/2

    list_of_points=[points]
    new_points=points

    if shape_of_matrix > 1:
        for i in range(1,int(shape_of_matrix)):
            n=2

            points_part1=new_points[0:n]
            points_part2=new_points[n:len(new_points)]

            changed_points=np.concatenate((points_part2, points_part1), axis=0)
            list_of_points.append(changed_points)
            new_points=changed_points
    else:
        list_of_points=[points]


    return list_of_points


def gauss_chebyshev(numerator_points, denumerator_points, limits, n=100):
    '''
    This function counts Gauss-Chebyshev integral
    '''
    x = np.cos((2*np.arange(n)+1)*np.pi/(2*n))*(limits[1]-limits[0])*0.5+np.mean(limits)
    y = np.ones(x.shape, np.complex)
    for p in numerator_points:
        y *= np.sqrt(np.abs(x-p))
        if limits[0]<=p:
            print ('rot -90')
            y *= 1j
    for p in denumerator_points:
        y /= np.sqrt(np.abs(x-p))
        if limits[0]<=p:
            print ('rot 90')
            y *= -1j
    return np.sum(y)*np.pi/n


class ConformalMapping:

    def __init__(self, points):

        self.points=np.asarray(points)

    def cl(self):
        '''
        Main part
        '''
        shape_of_matrix=(len(self.points)-2)//2
        print('Shape of matrix=', shape_of_matrix)

        #numerator_points, denumerator_points = create_numerator_and_denumerator_points(self.points)

        #list_of_points=function_for_points(self.points)
        list_=function_for_points(self.points)

        numerator_point_ids = [p for p in np.arange(2, len(self.points)-3, 2)]
        denumerator_point_ids = [p for p in np.arange(len(self.points)) if p not in numerator_point_ids]

        #print (numerator_point_ids, denumerator_point_ids)

        Q_mat = np.zeros(shape_of_matrix)
        phi_mat = np.zeros(shape_of_matrix)

        for i in range(int(shape_of_matrix)):
            integrals = []
            for j in range(len(self.points)-1):
                a_id, b_id = j, j+1

                numerator_gc = [p for p in numerator_point_ids]
                denumerator_gc = [p for p in denumerator_point_ids]

                if a_id in numerator_point_ids:
                    numerator_gc.append(a_id)
                else:
                    denumerator_gc.remove(a_id)

                if b_id in numerator_point_ids:
                    numerator_gc.append(b_id)
                else:
                    denumerator_gc.remove(b_id)


                numerator_gc_points = list_[i][numerator_gc]
                denumerator_gc_points = list_[i][denumerator_gc]
                limits = list_[i][[a_id, b_id]]

                integral = gauss_chebyshev(numerator_gc_points, denumerator_gc_points, limits)

                print (numerator_gc_points, denumerator_gc_points, limits)
                integrals.append(integral)
            transform = np.cumsum(integrals).tolist()
            print (transform)

            potentials = np.imag(transform[:-1:2]-trans)
            charges = np.real(integrals[1::2])

            #print (potentials, charges)
            #phi_mat[:, i] = potentials
            #Q_mat[:, i] = charges

        return list_of_points
