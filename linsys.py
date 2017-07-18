from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane
from para import Parameterization
from hyperplane import Hyperplane
getcontext().prec = 30


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'
    INF_SOLUTIONS = False

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
        self.planes[row] = self.planes[row].__class__(
            normal_vector=normal_vector, constant_term=constant_term)

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
        if indices[0] != -1:
            scalar = Decimal(
                1) / Decimal(rref.planes[0].normal_vector.coordinates[indices[0]])
            rref.multiply_coefficient_and_row(scalar, 0)
        return rref

    def GaussianEliminationSolution(self, esp=1e-10):
        solution = self.compute_rref()
        indices = solution.indices_of_first_nonzero_terms_in_each_row()
        pivot_variable_num = sum([1 if i >= 0 else 0 for i in indices])
        # for row in solution.planes[::-1]:
        # try:
        #     row.first_nonzero_index(row.normal_vector.coordinates)
        # except Exception as e:
        #     if str(e) == 'No nonzero elements found':
        #         constant_term = MyDecimal(row.constant_term)
        #         if not constant_term.is_near_zero():
        #             raise Exception(self.NO_SOLUTIONS_MSG)
        #     else:
        #         raise e

        for row in solution.planes[::-1]:
            if abs(sum(row.normal_vector.coordinates)) < esp and abs(row.constant_term) > esp:
                return self.NO_SOLUTIONS_MSG
            if pivot_variable_num < row.dimension:
                print(self.INF_SOLUTIONS_MSG)
                self.INF_SOLUTIONS = True
                para = self.InfiniteSolutionParamterization(
                    solution, row.dimension, pivot_variable_num, indices)
                return para
        # print(solution)
        return solution

    def InfiniteSolutionParamterization(self, rref, dimension, pivot_variable_num, indices):
        if not self.INF_SOLUTIONS:
            raise Exception('Not infinite solutions problem!')
        else:
            basepoint = [0] * dimension
            for index, coe in enumerate(indices):
                basepoint[coe] = rref.planes[index].constant_term
            direction_vector_list = list()
            for i in range(dimension):
                try:
                    indices.index(i)
                except Exception as e:
                    if str(e) == '{} is not in list'.format(i):
                        dir_vec = [0] * dimension
                        for row in range(min(len(indices), dimension)):
                            dir_vec[row] = - \
                                rref.planes[row].normal_vector.coordinates[i]
                        dir_vec[i] = 1
                        direction_vector_list.append(Vector(dir_vec))
            para = Parameterization(Vector(basepoint), direction_vector_list)
            return para


class MyDecimal(Decimal):

    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


if __name__ == '__main__':
    # p1 = Hyperplane(normal_vector=Vector([0.786, 0.786, 0.588]), constant_term=-0.714)
    # p2 = Hyperplane(normal_vector=Vector([-0.131, -0.131, 0.244]), constant_term=0.319)
    # s = LinearSystem([p1, p2])
    # r = s.GaussianEliminationSolution()
    # print(r)

    # p1 = Hyperplane(normal_vector=Vector([8.631, 5.112, -1.816]), constant_term=-5.113)
    # p2 = Hyperplane(normal_vector=Vector([4.315, 11.132, -5.27]), constant_term=-6.775)
    # p3 = Hyperplane(normal_vector=Vector([-2.158, 3.01, -1.727]), constant_term=-0.831)
    # s = LinearSystem([p1, p2, p3])
    # r = s.GaussianEliminationSolution()
    # print(r)

    # p1 = Hyperplane(normal_vector=Vector([0.935, 1.76, -9.365]), constant_term=-9.955)
    # p2 = Hyperplane(normal_vector=Vector([0.187, 0.352, -1.873]), constant_term=-1.991)
    # p3 = Hyperplane(normal_vector=Vector([0.374, 0.704, -3.746]), constant_term=-3.982)
    # p4 = Hyperplane(normal_vector=Vector([-0.561, -1.056, 5.619]), constant_term=5.973)
    # s = LinearSystem([p1, p2, p3, p4])
    # r = s.GaussianEliminationSolution()
    # print(r)

    # p1 = Hyperplane(normal_vector=Vector([0.786, 0.786, 8.123, 1.111, -8.363]), constant_term=-5.113)
    # p2 = Hyperplane(normal_vector=Vector([-0.131, 0.131, 7.05, -2.813, 1.19]), constant_term=-6.775)
    # p3 = Hyperplane(normal_vector=Vector([9.015, -5.873, -1.105, 2.013, -2.802]), constant_term=-0.831)
    # s = LinearSystem([p1, p2, p3])
    # r = s.GaussianEliminationSolution()
    # print(r)

    p1 = Hyperplane(normal_vector=Vector(
        [2.59892059451973e-10, 4.07250857161242e-07, 0.000638162093171666]), constant_term=53.2951176000000-1)
    p2 = Hyperplane(normal_vector=Vector(
        [2.07913647561578e-09, 1.62900342864497e-06, 0.00127632418634333]), constant_term=54.0022243800000-1)
    p3 = Hyperplane(normal_vector=Vector(
        [7.01708560520327e-09, 3.66525771445118e-06, 0.00191448627951500]), constant_term=54.7093311600000-1)
    p4 = Hyperplane(normal_vector=Vector(
        [1.66330918049263e-08, 6.51601371457987e-06, 0.00255264837268666]), constant_term=55.4164379400000-1)
    s = LinearSystem([p1, p2, p3, p4])
    r = s.GaussianEliminationSolution()
    print('Three Degree:')
    print(r)

    p1 = Hyperplane(normal_vector=Vector(
        [1.65853260658566e-13, 2.59892059451973e-10, 4.07250857161242e-07, 0.000638162093171666]), constant_term=53.2951176000000-1)
    p2 = Hyperplane(normal_vector=Vector(
        [2.65365217053706e-12, 2.07913647561578e-09, 1.62900342864497e-06, 0.00127632418634333]), constant_term=54.0022243800000-1)
    p3 = Hyperplane(normal_vector=Vector(
        [1.34341141133439e-11, 7.01708560520327e-09, 3.66525771445118e-06, 0.00191448627951500]), constant_term=54.7093311600000-1)
    p4 = Hyperplane(normal_vector=Vector(
        [4.24584347285929e-11, 1.66330918049263e-08, 6.51601371457987e-06, 0.00255264837268666]), constant_term=55.4164379400000-1)
    p5 = Hyperplane(normal_vector=Vector(
        [1.03658287911604e-10, 3.24865074314966e-08, 1.01812714290310e-05, 0.00319081046585833]), constant_term=56.1235447300000-1)
    s = LinearSystem([p1, p2, p3, p4, p5])
    r = s.GaussianEliminationSolution()
    print('Four Degree:')
    print(r)

    p1 = Hyperplane(normal_vector=Vector([1.05841263981216e-16, 1.65853260658566e-13, 2.59892059451973e-10, 4.07250857161242e-07, 0.000638162093171666]), constant_term=53.2951176000000-1)
    p2 = Hyperplane(normal_vector=Vector([3.38692044739892e-15, 2.65365217053706e-12, 2.07913647561578e-09, 1.62900342864497e-06, 0.00127632418634333]), constant_term=54.0022243800000-1)
    p3 = Hyperplane(normal_vector=Vector([2.57194271474356e-14, 1.34341141133439e-11, 7.01708560520327e-09, 3.66525771445118e-06, 0.00191448627951500]), constant_term=54.7093311600000-1)
    p4 = Hyperplane(normal_vector=Vector([1.08381454316766e-13, 4.24584347285929e-11, 1.66330918049263e-08, 6.51601371457987e-06, 0.00255264837268666]), constant_term=55.4164379400000-1)
    p5 = Hyperplane(normal_vector=Vector([3.30753949941301e-13, 1.03658287911604e-10, 3.24865074314966e-08, 1.01812714290310e-05, 0.00319081046585833]), constant_term=56.1235447300000-1)
    p6 = Hyperplane(normal_vector=Vector([8.23021668717939e-13, 2.14945825813502e-10, 5.61366848416262e-08, 1.46610308578047e-05, 0.00382897255902999]), constant_term=56.8306515100000-1)
    s = LinearSystem([p1, p2, p3, p4, p5, p6])
    r = s.GaussianEliminationSolution()
    print('Five Degree:')
    print(r)
