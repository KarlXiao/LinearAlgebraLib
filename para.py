from vector import Vector
from decimal import Decimal, getcontext

getcontext().prec = 30


class Parameterization(object):
    BASEPT_AND_DIR_VECTORS_MUST_BE_IN_THE_SAME_DIM_MSG = (
        'The basepoint and direction vectors should all live in the same dimensions!')

    def __init__(self, basepoint, direction_vectors):
        self.basepoint = basepoint
        self.direction_vectors = direction_vectors
        self.dimension = self.basepoint.dimension

        try:
            for v in direction_vectors:
                assert v.dimension == self.dimension
        except AssertionError:
            raise Exception(
                self.BASEPT_AND_DIR_VECTORS_MUST_BE_IN_THE_SAME_DIM_MSG)

    def __str__(self):
        def write_free_variables(index):
            output = ''
            if not index:
                output += 'x_{}'.format(index + 1)
            else:
                output += ', x_{}'.format(index + 1)
            return output

        def write_coordinates(coordinates):
            output = '['
            for i, j in enumerate(coordinates):
                if j >= 0:
                    output += '{}'.format(MyDecimal(abs(j)
                                                    ).zero_and_one_format())
                else:
                    output += '{}'.format(MyDecimal(-abs(j)
                                                    ).zero_and_one_format())
                output += ','
            return output[:-1] + ']'

        term = [write_free_variables(index) for index in range(self.dimension)]
        pivot = '[' + ''.join(term) + ']'

        dir_mat = ['+x_{}*'.format(x + 1) + write_coordinates(y.coordinates) for x, y in zip(
            self.indices_of_first_one_in_each_row(self.direction_vectors), self.direction_vectors)]
        dir_vec = ''.join(dir_mat)
        output = pivot + '=' + \
            write_coordinates(self.basepoint.coordinates) + dir_vec
        return output

    @staticmethod
    def indices_of_first_one_in_each_row(iterable):
        indices = [-1] * len(iterable)
        try:
            for ind, vec in enumerate(iterable):
                for k, c in enumerate(vec.coordinates):
                    if MyDecimal(c).is_near_one():
                        indices[ind] = k
        except Exception as e:
            raise e
        return indices


class MyDecimal(Decimal):

    def is_near_one(self, eps=1e-10):
        return abs(self - 1) < eps

    def zero_and_one_format(self, eps=1e-10):
        if abs(self) < eps:
            return 0
        elif abs(abs(self) - 1) < eps:
            return 1
        else:
            return self
if __name__ == "__main__":
    a = Vector([2, 1, 0, 0, 0])
    b = Vector([0, 0, 0, 1, 0])
    c = Vector([0, 0, 1, 0, 0])
    v1 = Parameterization(a, [b, c])
    print(v1)
