import numpy as np
import pylab
from matplotlib import pyplot as plot
from tabulate import tabulate

from ExecutionSummary import ExecutionSummary


def prepare_contour(levels=100, show_labels=False, fx=None):
    the_xs = []
    the_ys = []
    coord_minimum = -0.5
    coord_maximum = 1.3

    for val in np.arange(coord_minimum, coord_maximum*1.05 + 1e-100, 0.1):
        the_xs.append(val)
        the_ys.append(val)

    # Z values as a matrix
    the_values = np.ndarray((len(the_xs), len(the_ys)))

    # Populate values
    for the_x in range(0, len(the_xs)):
        for the_y in range(0, len(the_ys)):
            the_values[the_x][the_y] = fx(the_xs[the_x], the_ys[the_y])

    # axis limits
    pylab.xlim([coord_minimum, coord_maximum])
    pylab.ylim([coord_minimum, coord_maximum])

    plot.title('Contour plot')
    plot.xlabel('x')
    plot.ylabel('y')
    contours = plot.contour(the_xs, the_ys, the_values, levels=levels)

    # show values on contour
    if show_labels:
        plot.clabel(contours, inline=1, fontsize=10)
    plot.rcParams.update({'font.size': 12})


def summary_to_graph(s1, the_x_vector: [float, float], number_early_attempts=False, color="y-"):
    the_xs = []
    the_ys = []

    first_iteration = True
    i = -1
    for gamma, the_args, the_value in s1.gamma_x_value_history:
        i += 1
        if first_iteration:
            first_iteration = False

        the_xs.append(the_args[0])
        the_ys.append(the_args[1])

        if number_early_attempts:
            if 1 <= i+1 <= 4 and i > 0:
                plot.annotate(f'{i}', [the_xs[i], the_ys[i]])

    if number_early_attempts and [the_xs[1], the_ys[1]] != [the_x_vector[0], the_x_vector[1]]:
        plot.annotate(f'X', [the_x_vector[0], the_x_vector[1]])
    if len(the_xs) > 2 or not number_early_attempts:
        plot.annotate("*", [the_xs[-1], the_ys[-1]])

    plot.plot(the_xs, the_ys, color, label=s1.name)


def summary_simplex_to_graph(s1, number_early_attempts=False, color="r-", limit_steps=11111):
    vertices = []
    x_high_i = 0

    starting = True
    i = -1
    for gamma, the_args, the_value in s1.gamma_x_value_history:
        i += 1

        if starting:
            vertices.append(the_args)
            if i == 2:
                starting = False
            else:
                continue
        else:  # first occur with i == 3
            x_high_i = s1.simplex_high_history_indexes[i-3]
            vertices[x_high_i] = the_args

        # listened to change -> form triangle -> draw
        # form triangle
        triangle_xs = []
        triangle_ys = []
        for vertice in vertices:
            triangle_xs.append(vertice[0])
            triangle_ys.append(vertice[1])

        # triangle touches its start, doesn't it?
        triangle_xs.append(triangle_xs[0])
        triangle_ys.append(triangle_ys[0])

        # number
        if number_early_attempts:
            if 1 <= i+1 <= 7:
                if i == 2:
                    for iii in range(3):
                        plot.annotate(f'{iii+1}', [triangle_xs[iii], triangle_ys[iii]])
                else:
                    plot.annotate(f'{i+1}', [triangle_xs[x_high_i], triangle_ys[x_high_i]])

        # draw
        plot.plot(triangle_xs, triangle_ys, color, label=s1.name)
        if limit_steps == i:
            break
    plot.annotate("*", s1.gamma_x_value_history[-1][1])


def gamma_to_graph(the_summary: ExecutionSummary, color="r-"):
    the_ys = []
    for the_gamma, the_x, the_value in the_summary.gamma_x_value_history:
        the_ys.append(the_gamma)
    the_xs = range(len(the_ys))
    the_suffix = '.\nMetodas: ' + the_summary.name
    plot.title('Naudota gama reikšmė žingsnių juostoje' + the_suffix)
    plot.xlabel('žingsnis')
    plot.ylabel('gama')
    pylab.xlim(the_xs[0], the_xs[-1])
    the_max_y = max(*the_ys)
    the_max_y = the_max_y*1.1 if the_max_y*1.1 > 31 else 31
    pylab.ylim(0, the_max_y)
    plot.plot(the_xs, the_ys, color, label=the_summary.name + '. gamma')


def print_summary(*args: ExecutionSummary):
    headers = ["Starting point", "Start r", "End r",
               "Iterations",
               "Solution", "B(X, r)",
               "Steps (N)", "f calls", "df calls",
               "Volume (V)", "Backbone method"]

    if len(args) < 1:
        print(tabulate([], headers=headers))
        return
    # Use a breakpoint in the code line below to debug your script.
    data = []
    # print('Name, Solution, f(X), volume')
    for argument in args:
        data.append([argument.translated,
                     argument.r_start,
                     argument.r_end,
                     argument.r_iterations,
                     argument.solution,
                     argument.value,
                     argument.steps,
                     argument.fx_times,
                     argument.dfx_times,
                     argument.translated_fx,
                     argument.name]
                    )
        # print(f'{argument.name}, {argument.solution}, {argument.value}, {argument.translated_fx}')
    print('\n\n')
    print(tabulate(data, headers=headers))