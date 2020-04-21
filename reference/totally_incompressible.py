#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  9 15:41:43 2020.

@author: tagir.farkhutdinov@atco.com
"""

from collections import namedtuple
import warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import sympy as sy
from sympy.core import sympify
from sympy.utilities.autowrap import ufuncify


# %% Technical and term definitions
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


PROBLEM_CONSTS = {
    "rho_s": 1.,
    "rho_f": 1.,
    "g_0": 0.9,
    "K": 1,
    "param_a": 1.0,
    "param_b": 0,
    "S_0": 0.4,
    "nu": 0.0,
    "W": 1.,
    "U": 1.,
    "domain_half_length": 4.
}

FORMULAS = {
    "friction":
        sympify("K * (solid_t / density_s - fluid_t / density_f)"),


    "stress_tv":
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
                                     ),

    "stress_sponge":
        sympify('S_0 * exp( - (solid + x_coord)**2/W**2'
                ' - time/ U )'),

    "mu": (
        sympify("g_0*(X_t*Y_xt - X_x*R_Y) "
                "+ (1-g_0*Y_x)*((X_t*X_xt)/X_x - R_X)") /
        sympify("(X_x*Y_x*g_0)/rho_f + (1-g_0*Y_x)^2/rho_s")
        ),

    "mu_corr":
        sympify("X_x / "
                "((X_x*Y_x*g_0)/rho_f"
                "+ (g_0*Y_x - 1.)^2/rho_s)"),

    "p_x_times_g0_Y_x":
        sympify("param_b * g_0 * Y_x "
                " * (g_0*Y_xx + (1. - g_0)*X_xx)"),

    "sigma":
        sympify("- 0.5*param_a*(X_x**2 - 1.) "
                "- param_b*((X_x*(1. - g_0) + Y_x*g_0)**2 - 1.)"),

    "sigma_x":
        sympify('- param_a * X_x * X_xx'
                '- param_b'
                '* (g_0*Y_x + (1. - g_0)*X_x)'
                '* (g_0*Y_xx + (1. - g_0)*X_xx)'),

}


def ufunc_expr(symbols, formula, consts=PROBLEM_CONSTS):
    """Ufuncify the formula with symbols as arguments using numpy backend.

    Example: foo = ufunc_expr('x y', 'x % y')
    """
    if isinstance(symbols, str):
        symbols = sy.symbols(symbols)
    if isinstance(formula, str):
        formula = sympify(formula)
    if consts is not None:
        formula = formula.subs(
            [(sy.Symbol(k), v) for k, v in consts.items()]
            )
    return ufuncify(symbols, formula,
                    backend='numpy',
                    flags=['-D_USE_MATH_DEFINES'])


# Variable names
SYMB_ARGS = sy.symbols("time, x_coord, "
                       "solid, solid_t, density_s, laplace_s,"
                       "fluid, fluid_t, density_f, laplace_f")

# Alternative variable names
SY_ARG_XY = sy.symbols("t, x, "
                       "xi, X_t, X_x, X_xx,"
                       "eta, Y_t, Y_x, Y_xx")

# %% Interactive questions

stress = None

version = input("Please enter which stress to compile:\n"
                "'sponge' or 'tv'? >>>")
STRESS_TYPE = None
END_TIME = None
if version == 'tv':
    print("You've selected traveling wave stress")
    STRESS_TYPE = 'stress_tv'
    END_TIME = 50.
elif version == 'sponge':
    print("You've selected sponge stress")
    STRESS_TYPE = 'stress_sponge'
    END_TIME = 10.
else:
    print("Your selection was not recognized, assuming TV.")

end_time_new = input("Current END_TIME=%g. Enter new time\n"
                     "or keep current (press enter) >>>" % END_TIME)
if len(end_time_new):
    END_TIME = float(end_time_new)


RECOMPILE_UFUNCS = True
if RECOMPILE_UFUNCS:
    print("Compilation of ufuncs started...")
    global friction, pressure, sigma_x, sigma, mu, mu_corr, dynamic
    friction = ufunc_expr(SYMB_ARGS, FORMULAS['friction'])
    pressure = ufunc_expr(SY_ARG_XY,
                          FORMULAS['p_x_times_g0_Y_x'])
    sigma_x = ufunc_expr(SY_ARG_XY,
                         FORMULAS['sigma_x'])
    sigma = ufunc_expr(SY_ARG_XY,
                       FORMULAS['sigma'])

    stress = ufunc_expr(SYMB_ARGS, FORMULAS[STRESS_TYPE])

    mu = ufunc_expr('X_t, X_x, Y_x, X_xt, Y_xt, R_X, R_Y',
                    FORMULAS['mu'])

    mu_corr = ufunc_expr('X_t, X_x, Y_x, X_xt, Y_xt, R_X, R_Y',
                         FORMULAS['mu_corr'])

    dynamic = ufunc_expr('F_t, F_tx, density, laplace',
                         'F_t / density'
                         ' * (laplace * (F_t / density) - 2*F_tx)')
    print("Compilation of ufuncs succeeded!")
else:
    print("Ufuncs are already compiled!")


DIFF_SIGMA_TO_CRASH = False
if DIFF_SIGMA_TO_CRASH:
    warnings.warn("TF: You chose to differentiate Sigma numerically!",
                  UserWarning)


# %% Time step definition
def time_step(time, y, *args):
    """Compute one time step."""
    # pylint: disable=R0914
    x_coord = np.broadcast_to(args[0].reshape(-1, 1),
                              (args[0].size, y.shape[1]))

    solid, solid_t, fluid, fluid_t = np.split(y, 4)

    (rho_s, rho_f, g_0, K, param_a, param_b,
     S_0, nu, W, U, domain_half_length) = args[1]

    # domain_half_length = x_coord[-1, 0]
    partial_x = 2.*domain_half_length/(x_coord.size - 1)

    def apply_(arr, func, bc='s'):
        cols = range(arr.shape[1])
        return np.vstack([func(arr[:, i].T, bc) for i in cols]).T

    density_s = apply_(solid, partial, 's')/partial_x + 1.
    density_f = apply_(fluid, partial, 'f')/partial_x + 1.
    laplace_s = apply_(solid, partial2, 's')/partial_x**2
    laplace_f = apply_(fluid, partial2, 'f')/partial_x**2

    layer_data = (time, x_coord,
                  solid, solid_t, density_s, laplace_s,
                  fluid, fluid_t, density_f, laplace_f)

    problem_consts = ()

    friction_ = friction(*(layer_data + problem_consts))
    pressure_ = pressure(*(layer_data + problem_consts))

    sigma_ = None
    if not DIFF_SIGMA_TO_CRASH:
        sigma_ = sigma_x(*(layer_data + problem_consts))
    else:
        sigma_ = apply_(
            sigma(*(layer_data + problem_consts)),
            partial)/partial_x

    stress_ = stress(*(layer_data + problem_consts))
    stress_x = apply_(stress_, partial)/partial_x
    sigma_and_stress = stress_x + sigma_

    solid_xt = apply_(solid_t, partial, 's')/partial_x
    fluid_xt = apply_(fluid_t, partial, 'f')/partial_x

    delta_solid_t = (- dynamic(solid_t, solid_xt,
                               density_s, laplace_s)
                     - (friction_
                        + pressure_
                        + sigma_and_stress) / rho_s
                     )
    delta_fluid_t = (- dynamic(fluid_t, fluid_xt,
                               density_f, laplace_f)
                     + (friction_ + pressure_) / (g_0 * rho_f))

    incompr = mu(solid_t,
                 density_s, density_f,
                 solid_xt, fluid_xt,
                 delta_solid_t, delta_fluid_t)
    incompr_corr = mu_corr(solid_t,
                           density_s, density_f,
                           solid_xt, fluid_xt,
                           delta_solid_t, delta_fluid_t)
    incompr = incompr - incompr_corr * (
        incompr.mean(axis=0) / incompr_corr.mean(axis=0)
        )

    delta_solid_t += incompr * (density_f * -g_0 + 1.)/rho_s
    delta_fluid_t += incompr * density_f/rho_f

    return np.vstack((solid_t + nu * laplace_s,
                      delta_solid_t,
                      fluid_t + nu * laplace_f,
                      delta_fluid_t))

# %%


Statement = namedtuple('ProblemStatement',
                       ["number_of_intervals",
                        "time_interval",
                        "consts"])

Solution = namedtuple('ProblemSolution',
                      "x_coord, solid, solid_t, fluid, fluid_t")


def solve(problem_statement):
    """Solve PDE by reduction to IVP for ODE."""
    (number_of_intervals,
     time_interval,
     consts) = problem_statement
    domain_half_length = consts["domain_half_length"]
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
                       args=(x_coordinate,
                             np.array(list(consts.values())))
                       )

    if not isinstance(time_interval, np.ndarray):
        plt.plot(result.t)
        plt.title('Time (y) vs time steps (x).')
        plt.show()

    return Solution(x_coordinate,
                    *np.split(result.y.T, 4, axis=1))


def solve_instance(consts):
    """Solve problem instance."""
    # pylint: disable=W0612
    end_time = END_TIME
    statement = Statement(number_of_intervals=128,
                          time_interval=np.linspace(0., end_time,
                                                    int(end_time + 1.)),
                          consts=consts)

    return statement, solve(statement)


STATEMENT, SOLUTION = solve_instance(PROBLEM_CONSTS)

# %%


def plot_densities_and_velocities(solution, args={}):
    """Plot densities, physical and lagrangian velocities."""
    x_coord, solid, solid_t, fluid, fluid_t = solution
    partial_x = (x_coord[-1] - x_coord[0])/(x_coord.size - 1)
    time_index = solid.shape[0] - 1

    dens_s = partial(solid[time_index], 's')/partial_x + 1
    dens_f = partial(fluid[time_index], 'f')/partial_x + 1

    # plt.plot(x_coord, sigma_and_stress_term(0, solid + x_coordinate))
    plt.plot(x_coord, np.log(dens_s))
    plt.plot(x_coord, np.log(dens_f))
    plt.legend(('density of Solid',
                'density of Fluid'))
    plt.title('Logarithms of densities')
    plt.show()

    # Incomressibility layers, demonstrating C(t)
    for time_index in range(solid.shape[0] - 1):
        u_s = -solid_t[time_index] / (
            partial(solid[time_index], 's')/partial_x + 1.)
        u_f = -fluid_t[time_index] / (
            partial(fluid[time_index], 'f')/partial_x + 1.)

        dens_f = partial(fluid[time_index], 'f')/partial_x + 1

        g_0 = args['g_0']
        plt.plot(x_coord, g_0*dens_f*u_f + (1 - dens_f*g_0)*u_s)

    plt.legend((r'$gu_f+(1-g)u_s$', ))
    plt.title('Incomressibility expression')
    plt.show()

    plt.plot(x_coord, solid_t[time_index])
    plt.plot(x_coord, fluid_t[time_index])
    plt.legend((r'$X_t$',
                r'$Y_t$'))
    plt.title('Lagrangian Velocities')
    plt.show()


plot_densities_and_velocities(SOLUTION, STATEMENT.consts)

# %%


def plot_evolution(statement, solution):
    """Plot time dynamics."""
    number_of_intervals, time_interval, consts = statement
    domain_half_length = consts["domain_half_length"]
    x_coord, solid, solid_t, fluid, fluid_t = solution
    partial_x = (x_coord[-1] - x_coord[0])/(x_coord.size - 1)

    def plot2d():
        # plt.plot(x_coord, solid[-1])
        # plt.show()

        X, T = np.meshgrid(x_coord, time_interval)
        fig, (ax1, ax2) = plt.subplots(1, 2)
        # ax=fig.gca()
        # ax1=fig.gca(projection='3d')
        # ax.plot_surface(X,T,solid,cmap=cm.coolwarm)
        ax1.imshow(solid,
                   extent=[-domain_half_length, domain_half_length,
                           0, time_interval[-1]],
                   aspect='auto',
                   origin='lower')
        ax2.imshow(fluid,
                   extent=[-domain_half_length, domain_half_length,
                           0, time_interval[-1]],
                   aspect='auto',
                   origin='lower')
        # ax.view_init(90,0)
        ax1.set_xlabel('x')
        ax1.set_ylabel('t')
        ax1.set_title('Solid')
        ax2.set_xlabel('x')
        # ax2.set_ylabel('t')
        ax2.set_title('Fluid')

        # plt.savefig('plots/Evolution_XY.pdf')
        # plt.savefig('plots/Evolution_XY.png')
        plt.show()

    plot2d()

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

    rho_s, rho_f, g_0 = (statement.consts["rho_s"],
                         statement.consts["rho_f"],
                         statement.consts["g_0"],)

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

# %% Total incomressibility diagnostic test


def incompress_diag(solution, args={}):
    """Plot densities, physical and lagrangian velocities."""
    x_coord, solid, solid_t, fluid, fluid_t = solution
    partial_x = (x_coord[-1] - x_coord[0])/(x_coord.size - 1)
    time_index = solid.shape[0] - 1

    def apply_(arr, func):
        rows = range(arr.shape[0])
        return np.vstack([func(arr[i, :], 's') for i in rows])

    density_s = apply_(solid, partial)/partial_x + 1.
    density_f = apply_(fluid, partial)/partial_x + 1.

    plt.plot(x_coord,
             (1 - args['g_0'])*density_s[time_index]
             + args['g_0']*density_f[time_index])

    plt.legend('Sum of densities')
    plt.title('Sum of densities')
    plt.show()


incompress_diag(SOLUTION, STATEMENT.consts)
