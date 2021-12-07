from sympy import Symbol, lambdify
from sympy.parsing.sympy_parser import parse_expr

from Algorithms import *
from presentation import print_summary

if __name__ == '__main__':

    LSP = '1813056'
    LSP_a = int(LSP[4])
    LSP_b = int(LSP[5])
    LSP_c = int(LSP[6])

    variables = {'x': Symbol('x', real=True),
                 'y': Symbol('y', real=True),
                 'z': Symbol('z', real=True)
                 }
    variables_f = {
        'g1x': Symbol('g1x', real=True),
        'h1x': Symbol('h1x', real=True),
        'h2x': Symbol('h2x', real=True),
        'h3x': Symbol('h3x', real=True)
    }
    variables_bx = {
        'fx': Symbol('fx', real=True),
        'bx': Symbol('bx', real=True),
        'r': Symbol('r', real=True)
    }
    form_fx = parse_expr('-(x * y * z)', variables)
    form_g1x = parse_expr('2 * (x * y + x * z + y * z) - 1', variables)
    form_h1x = parse_expr('x * y * -1', variables)
    form_h2x = parse_expr('x * z * -1', variables)
    form_h3x = parse_expr('x * y * z * -1', variables)
    form_bx = parse_expr('g1x**2 + Max(0, h1x)**2 + Max(0, h2x)**2 + Max(0, h3x)**2', variables_f)
    form_penalty = parse_expr('fx + 1/r * bx', variables_bx)

    form_bx = form_bx.subs([
        (variables_f['g1x'], form_g1x),
        (variables_f['h1x'], form_h1x),
        (variables_f['h2x'], form_h2x),
        (variables_f['h3x'], form_h3x)
    ])
    form_penalty = form_penalty.subs([
        (variables_bx['bx'], form_bx),
        (variables_bx['fx'], form_fx)
    ])

    end_variable_values = [x for x in variables.values()] + [variables_bx['r']]

    fx = lambdify(variables.values(), form_fx, "numpy")
    Bx = lambdify(end_variable_values, form_penalty, "numpy")

    the_r = R()

    points = [[0, 0, 0],
              [1, 1, 1],
              [LSP_a / 10, LSP_b / 10, LSP_c / 10]
              ]

    epsilon = 1e-4
    the_r.value = 2
    the_r.multiplier = 0.2

    modules = [{'Heaviside': lambda x: np.heaviside(x, 1)}, 'numpy']
    form_dpenalty = [form_penalty.diff(the_symbol) for the_symbol in variables.values()]
    dpenalty = [lambdify(end_variable_values,
                         form_de_penalty,
                         modules)
                for form_de_penalty in form_dpenalty
                ]

    def gradient_penalty(values):
        return np.array([de_penalty(*values, the_r.value) for de_penalty in dpenalty])

    def call_penalty(*values):
        return Bx(*values, the_r.value)

    summaries = []
    for point in points:
        summaries.append(nonlinear_solve(point, the_r, call_penalty, gradient_penalty, fx, epsilon))

    print_summary(*summaries)

    # #    experiment. different r, r_mult, point
    # import warnings
    # warnings.filterwarnings("error")
    # summaries = []
    # r = [8, 6, 2, 1, 0.9, 0.8]
    # r_mult = [0.9, 0.8, 0.5, 0.2]
    # for rv in r:
    #     for rm in r_mult:
    #         the_r.value = rv
    #         the_r.multiplier = rm
    #         for point in points:
    #             try:
    #                 summaries.append(nonlinear_solve(point, the_r, call_penalty, gradient_penalty, fx, epsilon))
    #             except RuntimeWarning:
    #                 summaries.append(ExecutionSummary(name='?', solution='?', value='?', steps='?'))
    #             summaries[-1].r_start = str(summaries[-1].r_start) + f" rm: {rm}"
    # print_summary(*summaries)
