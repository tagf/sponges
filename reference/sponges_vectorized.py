#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 15:41:43 2020.

@author: tagir.farkhutdinov@atco.com
"""

from collections import namedtuple
import numpy as np
# import matplotlib
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import sympy as sy
from sympy.core import sympify
from sympy.utilities.autowrap import ufuncify


# %%
def partial(array, bc='s'):
    """Compute central discrete differences.

    The assumed inputs are X - x and Y - x.
    Example : partial X over x = partial(solid, 's') + 1.
    """
    result = np.correlate(array, np.array([-0.5, 0., 0.5]), mode='same')
    bc = 'periodic'
    if bc in ('f', 'fluid'):
        result[0] -= 0.5 * array[0]
        result[-1] += 0.5 * array[-1]
    elif bc in ('s', 'solid'):
        result[0] -= 0.
        result[-1] += 0.
    elif bc == 'periodic':
        result[0] -= 0.5 * array[-1]
        result[-1] += 0.5 * array[0]
    else:
        raise ValueError("Please specify valid boundary conditions.")
    return result


def partial2(array, bc='s'):
    """Compute the second discrete differences.
    The assumed inputs are X and Y.
    Example : partial2 X over x2 = partial2(solid, 's')/partial_x**2.
    """
    result = np.correlate(array, np.array([1., -2., 1.]), mode='same')
    bc = 'periodic'
    if bc in ('f', 'fluid'):
        result[0] += array[0]
        result[-1] += array[-1]
    elif bc in ('s', 'solid'):
        result[0] += 0.
        result[-1] += 0.
    elif bc == 'periodic':
        result[0] += array[-1]
        result[-1] += array[0]
    else:
        raise ValueError("Please specify valid boundary conditions.")
    return result


def ufunc_expr(symbols, formula):
    """Ufuncify the formula with symbols as arguments using numpy backend.

    Example: foo = ufunc_expr('x y', 'x % y')
    """
    if isinstance(symbols, str):
        symbols = sy.symbols(symbols)
    if isinstance(formula, str):
        formula = sympify(formula)
    return ufuncify(symbols, formula,
                    backend='numpy',
                    flags=['-D_USE_MATH_DEFINES'])


FORMULAS = {
    "friction":
        "K * (solid_t / density_s - fluid_t / density_f)",

    "pressure":  # TODO: Please confirm that the expression for the pressure term is correct (density_f)!
        "2 * param_b * g_0 * density_f "
        " * (g_0**2 * laplace_f + (1 - g_0) * laplace_s)",

    "sigma_x":
        '- 2. * param_a * density_s * laplace_s'
        '- 2. * param_b'
        '* (g_0**2*density_f + (1. - g_0)*density_s)'
        '* (g_0**2*laplace_f + (1. - g_0)*laplace_s)',


    "stress_v0":
        sympify('S_0 * exp( - X**2/W**2)'
                ).subs(sy.Symbol('X'),
                       sympify(
                           'sign(x)*(-domain_half_length + '
                           '(fabs(x) - domain_half_length) '
                           '% (2.*domain_half_length))'
                           )
                       ).subs(sympify('sign(x)**2'),
                              1
                              ).subs(sy.Symbol('x'),
                                     sympify(
                                         'solid + x_coord - U*time'
                                         )
                                     )

    }


# Variables
SYMB_ARGS = sy.symbols("time, x_coord, "
                       "solid, solid_t, density_s, laplace_s,"
                       "fluid, fluid_t, density_f, laplace_f")

# Constants
SYMB_ARGS += sy.symbols('rho_s, rho_f, g_0, K, '
                        'param_a, param_b, S_0, nu, W, U,'
                        'domain_half_length')

friction = ufunc_expr(SYMB_ARGS, FORMULAS['friction'])
pressure = ufunc_expr(SYMB_ARGS, FORMULAS['pressure'])
sigma_x = ufunc_expr(SYMB_ARGS, FORMULAS['sigma_x'])
stress = ufunc_expr(SYMB_ARGS, FORMULAS['stress_v0'])

dynamic = ufunc_expr('F_t, F_tx, density, laplace',
                     'F_t / density'
                     ' * (laplace * (F_t / density) - 2*F_tx)')


# %%

Statement = namedtuple('ProblemStatement',
                       ["number_of_intervals",
                        "domain_half_length",
                        "time_interval"])

Solution = namedtuple('ProblemSolution',
                      "x_coord, solid, solid_t, fluid, fluid_t")


def solve(problem_statement):
    """Solve PDE by reduction to IVP for ODE."""
    (number_of_intervals,
     domain_half_length,
     time_interval) = problem_statement
    x_coordinate = np.linspace(-domain_half_length,
                               domain_half_length,
                               number_of_intervals + 1)

    # Solid and fluid represent alpha and beta correspondingly
    solid = np.zeros_like(x_coordinate)
    fluid = np.zeros_like(x_coordinate)
    solid_t = np.zeros_like(solid)
    fluid_t = np.zeros_like(fluid)

    result = solve_ivp(time_step,
                       (time_interval[0], time_interval[-1]),
                       y0=np.hstack((solid, solid_t, fluid, fluid_t)),
                       method='LSODA',
                       t_eval=time_interval,  # time_interval,
                       vectorized=True,
                       args=(x_coordinate,)
                       )

    if not isinstance(time_interval, np.ndarray):
        plt.plot(result.t)
        plt.title('Time (y) vs time steps (x).')
        plt.show()

    return Solution(x_coordinate,
                    *np.split(result.y.T, 4, axis=1))


time_step_counter = 0
multipl_counter = 0


def time_step(time, y, *args):
    """Compute one time step."""
    # pylint: disable=R0914
    global time_step_counter, multipl_counter
    time_step_counter += 1
    multipl_counter += y.shape[1]

    solid, solid_t, fluid, fluid_t = np.split(y, 4)
    x_coord = np.broadcast_to(args[0].reshape(-1, 1),
                              (args[0].size, y.shape[1]))
    domain_half_length = x_coord[-1, 0]
    partial_x = 2.*domain_half_length/(x_coord.size - 1)
    rho_s = 1.
    rho_f = 1.
    g_0 = 0.5
    K = 1
    param_a = 1.
    param_b = 1.0
    S_0 = 0.1
    nu = 0.0
    W, U = 1, 1

    def apply_(arr, func):
        cols = range(arr.shape[1])
        return np.vstack([func(arr[:, i].T, 's') for i in cols]).T

    density_s = apply_(solid, partial)/partial_x + 1.
    density_f = apply_(fluid, partial)/partial_x + 1.
    laplace_s = apply_(solid, partial2)/partial_x**2
    laplace_f = apply_(fluid, partial2)/partial_x**2

    layer_data = (time, x_coord,
                  solid, solid_t, density_s, laplace_s,
                  fluid, fluid_t, density_f, laplace_f)

    problem_consts = (rho_s, rho_f, g_0, K,
                      param_a, param_b, S_0, nu, W, U,
                      domain_half_length)

    friction_ = friction(*(layer_data + problem_consts))
    pressure_ = pressure(*(layer_data + problem_consts))
    sigma_ = sigma_x(*(layer_data + problem_consts))

    stress_ = stress(*(layer_data + problem_consts))
    stress_x = apply_(stress_, partial)/partial_x
    sigma_and_stress = stress_x + sigma_

    delta_solid_t = (- dynamic(solid_t,
                               apply_(solid_t, partial)/partial_x,
                               density_s, laplace_s)
                     - (friction_
                        + pressure_
                        + sigma_and_stress) / rho_s
                     )
    delta_fluid_t = (- dynamic(fluid_t,
                               apply_(fluid_t, partial)/partial_x,
                               density_f, laplace_f)
                     + (friction_ + pressure_) / (g_0 * rho_f))

    return np.vstack((solid_t + nu * laplace_s,
                      delta_solid_t,
                      fluid_t + nu * laplace_f,
                      delta_fluid_t))


def solve_instance():
    """Solve problem instance."""
    # pylint: disable=W0612
    end_time = 480.
    statement = Statement(number_of_intervals=128,
                          domain_half_length=4.,
                          time_interval=np.linspace(0., end_time,
                                                    int(1.024 * end_time)))

    return statement, solve(statement)


STATEMENT, SOLUTION = solve_instance()

# %%


def plot_densities_and_velocities(solution):
    """Plot densities, physical and lagrangian velocities."""
    x_coord, solid, solid_t, fluid, fluid_t = solution
    partial_x = (x_coord[-1] - x_coord[0])/(x_coord.size - 1)
    time_index = solid.shape[0] - 1

    # plt.plot(x_coord, sigma_and_stress_term(0, solid + x_coordinate))
    plt.plot(x_coord, np.log(partial(solid[time_index], 's')/partial_x + 1))
    plt.plot(x_coord, np.log(partial(fluid[time_index], 'f')/partial_x + 1))
    plt.legend(('density of Solid',
                'density of Fluid'))
    plt.title('Logarithms of densities')
    plt.show()

    plt.plot(x_coord, -solid_t[time_index] /
             (partial(solid[time_index], 's')/partial_x + 1.))
    plt.plot(x_coord, -fluid_t[time_index] /
             (partial(fluid[time_index], 'f')/partial_x + 1.))
    plt.legend(('u_s',
                'u_f'))
    plt.title('Physical Velocities')
    plt.show()

    plt.plot(x_coord, solid_t[time_index])
    plt.plot(x_coord, fluid_t[time_index])
    plt.legend((r'$X_t$',
                r'$Y_t$'))
    plt.title('Lagrangian Velocities')
    plt.show()


plot_densities_and_velocities(SOLUTION)

# %%


def plot_evolution(statement, solution):
    """Plot time dynamics."""
    number_of_intervals, domain_half_length, time_interval = statement
    x_coord, solid, solid_t, fluid, fluid_t = solution
    partial_x = (x_coord[-1] - x_coord[0])/(x_coord.size - 1)
    if False:
        plt.plot(x_coord, solid[-1])
        plt.show()

        X, T = np.meshgrid(x_coord, time_interval)
        fig, (ax1, ax2) = plt.subplots(1, 2)
        # ax=fig.gca()
        # ax1=fig.gca(projection='3d')
        # ax.plot_surface(X,T,solid,cmap=cm.coolwarm)
        ax1.imshow(solid,
                   extent=[-domain_half_length, domain_half_length, 0, 10],
                   origin='lower')
        ax2.imshow(fluid,
                   extent=[-domain_half_length, domain_half_length, 0, 10],
                   origin='lower')
        # ax.view_init(90,0)
        ax1.set_xlabel('x')
        ax1.set_ylabel('t')
        ax1.set_title('Solid')
        ax2.set_xlabel('x')
        # ax2.set_ylabel('t')
        ax2.set_title('Fluid')

        # plt.show()
        plt.savefig('plots/Evolution_XY.pdf')
        plt.savefig('plots/Evolution_XY.png')

    # Traveling waves
    fig = plt.figure()
    ax = plt.axes()
    num_times = time_interval.size
    time_index_min = num_times - 10
    for j in range(time_index_min, num_times):
        time = time_interval[j]
        U = 1
        X1 = np.mod(x_coord-U*time, 2*domain_half_length) - domain_half_length
        ax.plot(X1, solid[j, :], 'o')
    ax.set_xlabel('x-U t')
    ax.set_ylabel('X(x,t)')
    plt.show()
    # plt.savefig('plots/Evolution.pdf')

    rho_s, rho_f, g_0 = (1, 1.0, 0.5)

    momentum_s = -rho_s*solid_t.sum(axis=1)*partial_x
    momentum_f = -g_0 * rho_f * fluid_t.sum(axis=1)*partial_x

    max_time_plot = num_times
    index_plot = range(max_time_plot)

    fig = plt.figure()
    ax = plt.axes()

    ax.plot(time_interval[index_plot],
            momentum_s[index_plot], 'b-')
    ax.plot(time_interval[index_plot],
            momentum_f[index_plot], 'r-')
    ax.plot(time_interval[index_plot],
            momentum_s[index_plot]+momentum_f[index_plot], 'k-')
    # plt.set_xlabel('t')
    ax.set_title('Momentum')
    ax.legend(('Solid', 'Fluid', 'Net'))
    plt.savefig('All_momenta.pdf')
    plt.savefig('All_momenta.png')
    ax.set_xlabel('t')
    ax.set_ylabel(r'$M_s$, $M_f$,$M_s+M_f$')

    plt.show()
    # fig.savefig('plots/Momentum.pdf')
    # fig.savefig('plots/Momentum.png')


plot_evolution(STATEMENT, SOLUTION)
