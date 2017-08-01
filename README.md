# LinearAlgebraLib

​	**Linear algebra library: solve linear system problems using gaussian elimination** 

##  `vector.py`

|           Function            | Description                              |
| :---------------------------: | :--------------------------------------- |
|           `__str__`           | print vector                             |
|           `__qr__`            | return true if $a$ equals $b$, false otherwise |
|           `__add__`           | return $a + b$                           |
|           `__sub__`           | return $a - b$                           |
|        `times_scalar`         | return $a*scalar$ ,same as `__mul__`     |
|          `magnitude`          | calculate the magnitude of vectors       |
|         `__truediv__`         | return $a/scalar$                        |
|         `unit_vector`         | calculate the unit vector                |
|         `dot_product`         | return dot product of $a$ and $b$        |
|            `angle`            | calculate the angle between `self` and another vector |
|           `is_zero`           | return true if it's a zero vector, false otherwise |
|      `is_orthogonal_to`       | return true if `self` is orthogonal to another vector, false otherwise |
|       `is_parallel_to`        | return true if `self` is paralle to another vector, false otherwise |
|        `projection_on`        | calculate the  prejection of `self` onto another vector |
|       ` orthogonal_on`        | calculate the component of `self` orthogonal to another vector |
|       ` cross_product`        | calculate the cross product of `self` and another vector |
| ` parallelogram_spanned_with` | calculate the area of parallelogram spanned by `self` and another vector |
|   ` triangle_spanned_with`    | calculate the area of triangle spanned by `self` and another vectro |

## `line.py`

|       Function        | Description                              |
| :-------------------: | ---------------------------------------- |
|       `__str__`       | print line                               |
| `first_nonzero_index` | find the first nonzero element of line normal vector |
|   ` is_paralle_to`    | return true if `self` is paralle to another line, false otherwise |
|       ` __eq__`       | return true if `self` equals another line, false otherwise |
| ` intersection_with`  | calculate the intersection of `self` and another |

## `plane.py`

|       Function        | Description                              |
| :-------------------: | ---------------------------------------- |
|       `__str__`       | print plane                              |
| `first_nonzero_index` | return the first nonzero element of plane normal vector |
|   ` is_paralle_to`    | return true if `self` ia paralle to another plane, false otherwise |
|       ` __eq__`       | return true if `self` equals to another plane, false otherwise |

## `hyperplane.py`

|       Function        | Description                              |
| :-------------------: | ---------------------------------------- |
|       `__str__`       | print hyperplane                         |
| `first_nonzero_index` | return the first nonzero element of hyperplane normal vector |
|   ` is_paralle_to`    | return true if `self` ia paralle to another hyperplane, false otherwise |
|       ` __eq__`       | return true if `self` equals to another hyperplane, false otherwise |

## `para.py`

​	**parameterize functions**

|              Function               | Description                              |
| :---------------------------------: | ---------------------------------------- |
|             ` __str__`              | print parameterized functions            |
| ` indices_of_first_one_in_each_row` | find every nonzero elements of direction vectors |

## `linsys.py`

​	**functions to operate the equation set(linear system)**

|                 Function                 | Description                              |
| :--------------------------------------: | ---------------------------------------- |
|               ` swap_rows`               | swap two rows in the equation set        |
|     ` multiply_coefficient_and_row`      | multiply a coefficient with indexed row  |
|     ` add_multiple_times_row_to_row`     | multiply a coefficient with indexed row, then plus another row |
| ` indices_of_first_nonzero_terms_in_each_row` | return indexes of first nonzeros elements of every row |
|                ` __len__`                | return total number of equations         |
|              ` __getitem__`              | return indexed equation                  |
|              ` __setitem__`              | set indexed equation to another          |
|                ` __str__`                | print equation set(linear system)        |
|        ` compute_triangular_form`        | compute triangular form of the equation set |
|             ` compute_rref`              | compute reduced row echelon form(rref) of the equation set |
|      ` GaussianEliminationSolution`      | apply gaussian elimination to rref of the equation set |
|    ` InfiniteSolutionParamterization`    | parameterize the gaussian elimination solution |