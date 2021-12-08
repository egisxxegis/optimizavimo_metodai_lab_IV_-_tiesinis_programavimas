from Algorithms import *

if __name__ == '__main__':

    LSP = '1813056'
    LSP_a = int(LSP[4])
    LSP_b = int(LSP[5])
    LSP_c = int(LSP[6])

    x_pratybos = \
        [[0,      3, -2, -6,  0, 0, 0, 0],
         [12,     4,  2,  0,  1, 1, 0, 0],
         [8,     -2,  2,  1, -1, 0, 1, 0],
         [5,      0,  0,  1,  2, 0, 0, 1]]

    x_uzduotis = \
        [[0,      2, -3,  0, -5, 0, 0, 0],
         [8,     -1,  1, -1, -1, 1, 0, 0],
         [10,     2,  4,  0,  0, 0, 1, 0],
         [3,      0,  0,  1,  1, 0, 0, 1]]

    x_mano_uzduotis = \
        [[0, 2, -3, 0, -5, 0, 0, 0],
         [LSP_a, -1, 1, -1, -1, 1, 0, 0],
         [LSP_b, 2, 4, 0, 0, 0, 1, 0],
         [LSP_c, 0, 0, 1, 1, 0, 0, 1]]

    solve_linear(x_uzduotis)
    solve_linear(x_mano_uzduotis)
    # solve_linear(x1)
