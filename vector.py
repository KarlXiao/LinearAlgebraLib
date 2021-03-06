import math
from decimal import Decimal, getcontext

getcontext().prec = 30

class Vector(object):

    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple([Decimal(x) for x in coordinates])
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')

    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)

    def __eq__(self, v):
        if self and v:
            return self.coordinates == v.coordinates
        else:
            return self.coordinates == None and v.coordinates == None

    def __add__(self, v):
        ''' calculate vector a plus vector b'''
        if self.dimension == v.dimension:
            out = [x + y for x, y in zip(self.coordinates, v.coordinates)]
            return self.__class__(out)
        else:
            print('Two vectors must have the same dimensions!')

    def __sub__(self, v):
        ''' calculate vector a minus vector b'''
        if self.dimension == v.dimension:
            out = [x - y for x, y in zip(self.coordinates, v.coordinates)]
            return self.__class__(out)
        else:
            print('Two vectors must have the same dimensions!')

    def __mul__(self, scalar):
        '''scalling vectors'''
        out = [x * Decimal(scalar) for x in self.coordinates]
        return self.__class__(out)

    def __rmul__(self, scalar):
        '''scalling vectors'''
        out = [x * Decimal(scalar) for x in self.coordinates]
        return self.__class__(out)

    def times_scalar(self,scalar):
        '''scalling vectors'''
        out = [x * Decimal(scalar) for x in self.coordinates]
        return self.__class__(out)

    def magnitude(self):
        '''calculate the magnitude of vectors'''
        out = Decimal(math.sqrt(sum(x**2 for x in self.coordinates)))
        return out

    def __truediv__(self, scalar):
        try:
            out = [x / Decimal(scalar) for x in self.coordinates]
            return self.__class__(out)
        except ZeroDivisionError:
            raise Exception('Cannot divided by zero!')

    def unit_vector(self):
        '''calculate the unit vector'''
        return self / self.magnitude()

    def dot_product(self, v):
        '''calculate the dot product of vector self and vector v'''
        if self.dimension == v.dimension:
            out = sum([x * y for x, y in zip(self.coordinates, v.coordinates)])
            return out
        else:
            print('Tow vectors must have the same dimensions!')

    def angle(self, v, mode='rad'):
        '''calculate the angle between vector self and vector v'''
        try:
            out = Decimal(math.acos(max(min(self.dot_product(
                v) / (self.magnitude() * v.magnitude()), 1), -1)) % (2 * math.pi))
            if mode == 'rad':
                return out
            else:
                return out / Decimal(math.pi * 180)
        except ZeroDivisionError:
            raise Exception('Vectors cannot be zero!')

    def is_zero(self, tolerance=1e-10):
        '''judge the vector is 0 or not'''
        return abs(self.magnitude()) < tolerance

    def is_orthogonal_to(self, v, tolerance=1e-10):
        '''judge vector self is orthogonal to vector v or not'''
        return abs(self.dot_product(v)) < tolerance

    def is_parallel_to(self, v, tolerance=1e-6):
        '''judge vector self is parallel to vector v or not'''
        return self.is_zero() or v.is_zero() or abs(self.angle(v) - Decimal(math.pi)) < tolerance or abs(self.angle(v)) < tolerance

    def projection_on(self, v):
        '''calculate the prejection of self onto v'''
        unit_vector = v.unit_vector()
        return unit_vector * self.dot_product(unit_vector)

    def orthogonal_on(self, v):
        '''calculate the component of self orthogonal to v'''
        return self - self.projection_on(v)

    def cross_product(self, v):
        '''calculate the cross product of vector self and vector v'''
        try:
            x_1, y_1, z_1 = self.coordinates
            x_2, y_2, z_2 = v.coordinates
            return self.__class__([y_1 * z_2 - y_2 * z_1, -(x_1 * z_2 - x_2 * z_1), x_1 * y_2 - x_2 * y_1])
        except ValueError as e:
            msg = str(e)
            if msg == 'need more than 2 values to unpack':
                self_embedded_in_threeDimension = Vector(
                    self.coordinates + ('0',))
                v_embedded_in_threeDimension = Vector(v.coordinates + ('0',))
                return self_embedded_in_threeDimension.cross_product(v_embedded_in_threeDimension)
            elif(msg == 'too many values to unpack' or msg == 'need more than 1 value to unpack'):
                raise Exception("Cross product only defined in tow three dimensions vectors!")
            else:
                raise e

    def parallelogram_spanned_with(self, v):
        '''calculate the area of parallelogram spanned by vector self and vector v'''
        # return self.cross_product(v).magnitude()
        return self.magnitude() * v.magnitude() * Decimal(math.sin(self.angle(v)))

    def triangle_spanned_with(self, v):
        '''calculate the area of triangle spanned by vector self and vectro v'''
        return Decimal(0.5) * self.parallelogram_spanned_with(v)

if __name__ == '__main__':
    v1 = Vector([8.218, -9.341])
    v2 = Vector([-1.129, 2.111])
    print('Add:', v1 + v2)

    v3 = Vector([7.119, 8.215])
    v4 = Vector([-8.223, 0.878])
    print('Minus:', v3 - v4)

    v5 = Vector([1.671, -1.012, -.318])
    a = 7.41
    print('Scalling:', v5 * a)

    v6 = Vector([-0.221, 7.437])
    print('Magnitude:', v6.magnitude())

    v7 = Vector([5.581, -2.136])
    print('Unit Vector:', v7.unit_vector())

    v8 = Vector([-5.955, -4.904, -1.874])
    v9 = Vector([-4.496, -8.755, 7.103])
    print(v8.dot_product(v9))

    v10 = Vector([7.35, 0.221, 5.188])
    v11 = Vector([2.751, 8.259, 3.985])
    print(v10.angle(v11, 'degree'))

    v12 = Vector([-2.328, -7.284, -1.214])
    v13 = Vector([-1.821, 1.072, -2.94])
    print(v12.is_parallel_to(v13))
    print(v12.is_orthogonal_to(v13))

    v14 = Vector([3.009, -6.172, 3.692, -2.51])
    v15 = Vector([6.404, -9.144, 2.759, 8.718])
    print(v14.projection_on(v15))
    print(v14.orthogonal_on(v15))

    v16 = Vector([8.462, 7.893, -8.187])
    v17 = Vector([6.984, -5.975, 4.778])
    print(v16.cross_product(v17))

    v16 = Vector([-8.987, -9.838, 5.031])
    v17 = Vector([-4.268, -1.861, -8.866])
    print(v16.parallelogram_spanned_with(v17))

    v16 = Vector([1.5, 9.547, 3.691])
    v17 = Vector([-6.007, 0.124, 5.772])
    print(v16.triangle_spanned_with(v17))

    v18 = Vector([1.671, -1.012, -.318])
    a = 7.41
    print('Scalling:', a * v18)