"""Modifications of potential."""

def pde_info(_sp=sympy):
    """ insert the below code into the function """

    print('Potential energy:')
    rho_s, rho_f = X.diff(x), g/g_0

    """Vakhtang's new potential:"""
    def potential_rhs(rho_s, rho_f):
        old = _sp.Symbol('alpha')/2 * (rho_s - 1)**2 + _sp.Symbol('beta')/2 * (rho_f*g_0 - (1-(1-g_0)*rho_s))**2
        return  _sp.Symbol('alpha')/2 * (1 - g_0 * rho_f) * (rho_s-1)**2 

    eq_potential = _sp.Eq(_sp.Symbol('V'),
                          potential_rhs(rho_s, rho_f)
                          )
    print_eqn(eq_potential, True)
    
