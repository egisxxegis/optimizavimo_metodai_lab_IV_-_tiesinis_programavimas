import copy
import numpy as np


def init_matrix(x: [[float]]):
    if len(x) < 1:
        x = [[0, 3, -2, -6, 0,	0, 0, 0],
             [12, 4, 2, 0, 1, 1, 0, 0],
             [8, -2, 2, 1, -1, 0, 1, 0],
             [5, 0, 0, 1, 2, 0, 0, 1]]
    return np.array(x, dtype='float64')


def matrix_line_sum_to_destination(matrix: np.ndarray, line_number_from_top_source: int, origin_multiplier: float,
                                   line_number_from_top_destination: int):
    if not isinstance(matrix, np.ndarray):
        raise TypeError("Duota matrix nebuvo numpy.ndarray tipo")

    source = matrix[line_number_from_top_source]
    destination = matrix[line_number_from_top_destination]
    matrix[line_number_from_top_destination] = source * origin_multiplier + destination
    return matrix


def matrix_line_multiplication(matrix: np.ndarray, line_number_from_top: int, multiplier: float):
    if not isinstance(matrix, np.ndarray):
        raise TypeError("Duota matrix nebuvo numpy.ndarray tipo")

    matrix[line_number_from_top] = matrix[line_number_from_top] * multiplier
    return matrix


def extract_col(matrix: np.ndarray, col_num, ignore_line_count=0):
    if not isinstance(matrix, np.ndarray):
        raise TypeError("Duota matrix nebuvo numpy.ndarray tipo")

    return matrix[ignore_line_count:, col_num]


def cols_to_matrix(*args):
    return np.stack(args, axis=1)


def find_bases_ind(matrix: np.ndarray, count_of_top_rows_to_be_ignored=1, cols0=None):

    cols_ind = [x for x in range(1, len(matrix)+1 - count_of_top_rows_to_be_ignored)] if cols0 is None else cols0

    def current_det():
        copied = copy.deepcopy(cols_ind)
        copied.sort()
        cols = [extract_col(matrix, x, count_of_top_rows_to_be_ignored) for x in copied]
        return np.linalg.det(cols_to_matrix(*cols))

    if current_det() >= 0:
        return cols_ind

    ind_max = len(matrix[0]) - 1
    for lead_i in range(0, ind_max + 1 - len(cols_ind)):
        for i in range(1, len(cols_ind)):   # for each not first index of indexes
            start = cols_ind[i]
            for ii in range(start, ind_max+1):   # try each following index
                if cols_ind.count(ii) > 0:
                    continue
                cols_ind[i] = ii
                if current_det() >= 0:
                    return cols_ind.sort()
            cols_ind[i] = start
        cols_ind[0] = cols_ind[-1] + 1
        if current_det() >= 0:
            return cols_ind.sort()
        cols_ind.sort()
        # continue search
    # after for loop
    raise TypeError("Duota matrix neturi bazės")


def find_base_row(matrix: np.ndarray, col_ind, count_of_top_rows_to_be_ignored=1):
    col = extract_col(matrix, col_ind, count_of_top_rows_to_be_ignored)
    for i in range(len(col)):
        if col[i] == 1:
            return i + count_of_top_rows_to_be_ignored
    return None


class Base:
    def __init__(self, col, row):
        self.col = col
        self.row = row


def solve_linear(matrix0: [[float]], ignore_the_top_lines_count=1):
    ignore_top = ignore_the_top_lines_count
    matrix = init_matrix(copy.deepcopy(matrix0))
    temp_bases = find_bases_ind(matrix, count_of_top_rows_to_be_ignored=ignore_top,
                                cols0=[-3 + x for x in range(len(matrix) - 1)])  # [-3-2-1]
    bases = [Base(x, find_base_row(matrix, x, ignore_top)) for x in temp_bases]
    for base in bases:
        translated_col = base.col if base.col > -1 else len(matrix[0]) + base.col
        base.col = translated_col

    print("Turime tokią matricą")
    print(matrix)

    def find_lowest_negative_top_ind(start_ind=1):
        the_min = matrix[0][start_ind]
        the_ind = start_ind
        for i in range(start_ind, len(matrix[0])):
            if matrix[0][i] < the_min:
                the_ind = i
                the_min = matrix[0][i]
        return the_ind if the_min < 0 else None

    def calc_lambdas_for_change(col_ind, ignore_top_lines=ignore_top):
        left_values = extract_col(matrix, 0, ignore_top_lines)
        right_values = extract_col(matrix, col_ind, ignore_top_lines)
        right_values_mask = right_values > 0
        right_values = right_values * right_values_mask
        right_values[right_values == 0] = 1e-12
        left_values = left_values / right_values
        left_values[left_values == 0] = np.Infinity
        return left_values

    def pick_lowest_ind(col_matrix, ignore_top_lines_from_earlier=ignore_top):
        the_min = col_matrix[0]
        the_ind = 0
        for i in range(len(col_matrix)):
            if col_matrix[i] < the_min:
                the_min = col_matrix[i]
                the_ind = i
        return the_ind + ignore_top_lines_from_earlier

    def transfer_into_base(col_ind, row_ind):
        the_new_base = Base(col_ind, row_ind)
        for i in range(len(bases)):
            if (bases[i]).row == the_new_base.row:
                bases[i] = the_new_base
                return the_new_base
        return None

    def recalculate_matrix_by_new_base(_base: Base):
        # normalise row
        base_value = matrix[_base.row][_base.col]
        multiplier = 1 / base_value
        if multiplier != 1:
            print(f"Dauginame eilutę {_base.row+1} (nuo apačios eilutę {len(matrix)-_base.row}) "
                  f"iš {multiplier}; kitaip:")
            print(f"Eilutė{_base.row + 1} = Eilutė{_base.row + 1} * {multiplier}\n")
            matrix_line_multiplication(matrix, _base.row, multiplier)
            base_value = matrix[_base.row][_base.col]
            print("Gauname:")
            print(matrix)

        # zeroing other rows
        for i in range(len(matrix)):
            if i == _base.row:
                continue
            if matrix[i][_base.col] == 0:
                continue
            multiplier = matrix[i][_base.col] / base_value * -1
            print(f"Padauginsime eilutę {_base.row + 1} (nuo apačios eilutę {len(matrix) - _base.row}) "
                  f"iš {multiplier} ir"
                  f" pridedame prie eilutės {i+1}; kitaip:")
            print(f"Eilutė{i+1} = Eilutė{i+1} + Eilutė{_base.row + 1} * {multiplier}\n")
            matrix_line_sum_to_destination(matrix, _base.row, multiplier, i)
            print("Gauname:")
            print(matrix)

    # main
    iteration = 0
    while True:
        iteration += 1
        new_base_ind = find_lowest_negative_top_ind(1)
        if new_base_ind is None:
            values = extract_col(matrix, 0, 0)
            variables = [0] * (len(matrix[0]) - len(matrix))
            for base in bases:
                translated_col = base.col if base.col > -1 else len(matrix[0]) + base.col
                if translated_col - 1 < len(variables):
                    variables[translated_col-1] = values[base.row]
            bases.sort(key=lambda x: x.col)
            print(f"Radome sprendinį X = {variables} su reikšme {values[0] * -1}")
            print(f"Bazių indeksai: {[x.col for x in bases]}")
            print("-----------------------------------\n")
            return
        lambdas = calc_lambdas_for_change(new_base_ind)
        new_base_row_ind = pick_lowest_ind(lambdas)
        new_base = transfer_into_base(new_base_ind, new_base_row_ind)
        recalculate_matrix_by_new_base(new_base)
        print(f"--------------Iteracija {iteration} baigta")
