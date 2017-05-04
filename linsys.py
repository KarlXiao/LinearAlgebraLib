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
        self[row1], self[row2] = self[row2], self[row1]

    def multiply_coefficient_and_row(self, coefficient, row):
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
                indices[i] = p.first_nonzero_index(p.normal_vector.coordinates)
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
        for row_idx in range(len(system) - 1):
            indices = system.indices_of_first_nonzero_terms_in_each_row()
            first_nonzero_row = row_idx + 1
            while row_idx != indices[row_idx] and first_nonzero_row < len(system):
                if row_idx == indices[first_nonzero_row]:
                    system.swap_rows(row_idx, first_nonzero_row)
                    indices = system.indices_of_first_nonzero_terms_in_each_row()
                first_nonzero_row += 1
            if row_idx == indices[row_idx]:
                for row_to_compute in range(row_idx + 1, len(system)):
                    coefficient = system.planes[row_to_compute].normal_vector.coordinates[
                        row_idx] / system.planes[row_idx].normal_vector.coordinates[row_idx] * Decimal(-1)
                    system.add_multiple_times_row_to_row(
                        coefficient, row_idx, row_to_compute)
        return system

    def compute_rref(self):
        rref = self.compute_triangular_form()
        for row_idx in range(len(rref))[::-1]:
            indices = rref.indices_of_first_nonzero_terms_in_each_row()
            if indices[row_idx] != -1:
                for row_to_compute in range(row_idx):
                    coefficient = rref.planes[row_to_compute].normal_vector.coordinates[indices[
                        row_idx]] / rref.planes[row_idx].normal_vector.coordinates[indices[row_idx]] * Decimal(-1)
                    rref.add_multiple_times_row_to_row(
                        coefficient, row_idx, row_to_compute)
                    scalar = Decimal(
                        1) / Decimal(rref.planes[row_idx].normal_vector.coordinates[indices[row_idx]])
                    rref.multiply_coefficient_and_row(scalar, row_idx)
        return rref

    def GaussianEliminationSolution(self):
        solution = self.compute_rref()
        indices = solution.indices_of_first_nonzero_terms_in_each_row()
        for row in range(len(solution))[::-1]:
            if sum(solution.planes[row].normal_vector.coordinates) == 0 and solution.planes[row].constant_term !=0:
                print('No Solution!')
            elif:
                pass
                print(solution + '\n')
                print('Infinite Solutions!')
            else:
                print(solution)
        return solution

class MyDecimal(Decimal):

    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps

if __name__ == '__main__':
    p1 = Plane(normal_vector=Vector([1, 1, 1]), constant_term=1)
    p2 = Plane(normal_vector=Vector([0, 1, 1]), constant_term=2)
    s = LinearSystem([p1, p2])
    r = s.compute_rref()
    print(r)
    if not (r[0] == Plane(normal_vector=Vector([1, 0, 0]), constant_term=-1) and
            r[1] == p2):
        print('test case 1 failed')

    p1 = Plane(normal_vector=Vector([1, 1, 1]), constant_term=1)
    p2 = Plane(normal_vector=Vector([1, 1, 1]), constant_term=2)
    s = LinearSystem([p1, p2])
    r = s.compute_rref()
    print(r)
    if not (r[0] == p1 and
            r[1] == Plane(constant_term=1)):
        print('test case 2 failed')

    p1 = Plane(normal_vector=Vector([1, 1, 1]), constant_term=1)
    p2 = Plane(normal_vector=Vector([0, 1, 0]), constant_term=2)
    p3 = Plane(normal_vector=Vector([1, 1, -1]), constant_term=3)
    p4 = Plane(normal_vector=Vector([1, 0, -2]), constant_term=2)
    s = LinearSystem([p1, p2, p3, p4])
    r = s.compute_rref()
    print(r)
    if not (r[0] == Plane(normal_vector=Vector([1, 0, 0]), constant_term=0) and
            r[1] == p2 and
            r[2] == Plane(normal_vector=Vector([0, 0, -2]), constant_term=2) and
            r[3] == Plane()):
        print('test case 3 failed')

    p1 = Plane(normal_vector=Vector([0, 1, 1]), constant_term=1)
    p2 = Plane(normal_vector=Vector([1, -1, 1]), constant_term=2)
    p3 = Plane(normal_vector=Vector([1, 2, -5]), constant_term=3)
    s = LinearSystem([p1, p2, p3])
    r = s.compute_rref()
    print(r)
    if not (r[0] == Plane(normal_vector=Vector([1, 0, 0]), constant_term=Decimal(23) / Decimal(9)) and
            r[1] == Plane(normal_vector=Vector([0, 1, 0]), constant_term=Decimal(7) / Decimal(9)) and
            r[2] == Plane(normal_vector=Vector([0, 0, 1]), constant_term=Decimal(2) / Decimal(9))):
        print('test case 4 failed')
