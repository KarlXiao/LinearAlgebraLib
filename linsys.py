from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane

getcontext().prec = 30


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    def swap_rows(self, row1, row2):
        # plane1 = self.planes[row1]
        # self.planes[row1] = self.planes[row2]
        # self.planes[row2] = plane1
        self[row1], self[row2] = self[row2], self[row1]

    def multiply_coefficient_and_row(self, coefficient, row):
        # self.planes[row].constant_term = self.planes[
        #     row].constant_term * coefficient
        # self.planes[row].normal_vector = [
        #     x * coefficient for x in self.planes[row].normal_vector]
        constant_term = self.planes[row].constant_term * coefficient
        normal_vector = self.planes[row].normal_vector * coefficient
        self.planes[row] = Plane(normal_vector, constant_term)

    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        self.planes[
            row_to_be_added_to].constant_term += self.planes[row_to_add].constant_term * coefficient
        self.planes[row_to_be_added_to].normal_vector = self.planes[
            row_to_add].normal_vector * coefficient + self.planes[row_to_be_added_to].normal_vector
        self.planes[row_to_be_added_to].set_basepoint()

    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i, p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices

    def __len__(self):
        return len(self.planes)

    def __getitem__(self, i):
        return self.planes[i]

    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i + 1, p)
                for i, p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret

    def compute_triangular_form(self):
        system = deepcopy(self)


class MyDecimal(Decimal):

    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps

if __name__ == '__main__':
    p1 = Plane(Vector([1, 1, 1]), 1)
    p2 = Plane(Vector([0, 1, 1]), 2)
    s = LinearSystem([p1, p2])
    t = s.compute_triangular_form()
    if not (t[0] == p1 and
            t[1] == p2):
        print('test case 1 failed')

    p1 = Plane(Vector([1, 1, 1]), 1)
    p2 = Plane(Vector([1, 1, 1]), 2)
    s = LinearSystem([p1, p2])
    t = s.compute_triangular_form()
    if not (t[0] == p1 and
            t[1] == Plane(1)):
        print('test case 2 failed')

    p1 = Plane(Vector([1, 1, 1]), 1)
    p2 = Plane(Vector([0, 1, 0]), 2)
    p3 = Plane(Vector([1, 1, -1]), 3)
    p4 = Plane(Vector([1, 0, -2]), 2)
    s = LinearSystem([p1, p2, p3, p4])
    t = s.compute_triangular_form()
    if not (t[0] == p1 and
            t[1] == p2 and
            t[2] == Plane(Vector([0, 0, -2]), 2) and
            t[3] == Plane()):
        print('test case 3 failed')

    p1 = Plane(Vector([0, 1, 1]), 1)
    p2 = Plane(Vector([1, -1, 1]), 2)
    p3 = Plane(Vector([1, 2, -5]), 3)
    s = LinearSystem([p1, p2, p3])
    t = s.compute_triangular_form()
    if not (t[0] == Plane(Vector([1, -1, 1]), 2) and
            t[1] == Plane(Vector([0, 1, 1]), 1) and
            t[2] == Plane(Vector([0, 0, -9]), -2)):
        print('test case 4 failed')
