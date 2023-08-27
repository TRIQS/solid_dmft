#!@TRIQS_PYTHON_EXECUTABLE@

import numpy as np
import collections
from itertools import product

'''
python functions for reading Uijkl from a VASP cRPA run and the evaluating the matrix
elements for different basis sets.

Copyright (C) 2020, A. Hampel and M. Merkel from Materials Theory Group
at ETH Zurich
'''


def read_uijkl(path_to_uijkl, n_sites, n_orb):
    '''
    reads the VASP UIJKL files or the vijkl file if wanted

    Parameters
    ----------
    path_to_uijkl : string
        path to Uijkl like file
    n_sites: int
        number of different atoms (Wannier centers)
    n_orb : int
        number of orbitals per atom

    Returns
    -------
    uijkl : numpy array
        uijkl Coulomb tensor

    '''
    dim = n_sites*n_orb
    uijkl = np.zeros((dim, dim, dim, dim), dtype=complex)
    data = np.loadtxt(path_to_uijkl)

    for line in range(0, len(data[:, 0])):
        i = int(data[line, 0])-1
        j = int(data[line, 1])-1
        k = int(data[line, 2])-1
        l = int(data[line, 3])-1
        uijkl[i, j, k, l] = data[line, 4]+1j*data[line, 5]

    return uijkl


def construct_U_kan(n_orb, U, J, Up=None, Jc=None):
    '''
    construct Kanamori Uijkl tensor for given U, J, Up, and Jc

    Parameters
    ----------
    n_orb : int
        number of orbitals
    U : float
        U value for elements Uiiii
    J : float
        Hunds coupling J for tensor elements Uijji
    Up : float, optional, default=U-2J
        inter orbital exchange term Uijij
    Jc : float, optional, default=J
        Uiijj term, is the same as J for real valued wave functions

    Returns
    -------
    uijkl : numpy array
        uijkl Coulomb tensor

    '''

    orb_range = range(0, n_orb)
    U_kan = np.zeros((n_orb, n_orb, n_orb, n_orb))

    if not Up:
        Up = U-2*J
    if not Jc:
        Jc = J

    for i, j, k, l in product(orb_range, orb_range, orb_range, orb_range):
        if i == j == k == l:  # Uiiii
            U_kan[i, j, k, l] = U
        elif i == k and j == l:  # Uijij
            U_kan[i, j, k, l] = Up
        elif i == l and j == k:  # Uijji
            U_kan[i, j, k, l] = J
        elif i == j and k == l:  # Uiijj
            U_kan[i, j, k, l] = Jc
    return U_kan


def red_to_2ind(uijkl, n_sites, n_orb, out=False):
    '''
    reduces the 4index coulomb matrix to a 2index matrix and
    follows the procedure given in PRB96 seth,peil,georges:
    U_antipar = U_mm'^oo' = U_mm'mm' (Coulomb Int)
    U_par = U_mm'^oo = U_mm'mm' - U_mm'm'm (for intersite interaction)
    U_ijij (Hunds coupling)
    the indices in VASP are switched: U_ijkl ---VASP--> U_ikjl

    Parameters
    ----------
    uijkl : numpy array
        4d numpy array of Coulomb tensor
    n_sites: int
        number of different atoms (Wannier centers)
    n_orb : int
        number of orbitals per atom
    out : bool
        verbose mode

    Returns
    -------
    Uij_anti : numpy array
        red 2 index matrix U_mm'mm'
    Uiijj : numpy array
        red 2 index matrix U_iijj
    Uijji : numpy array
        red 2 index matrix Uijji
    Uij_par : numpy array
        red 2 index matrix U_mm\'mm\' - U_mm\'m\'m
    '''
    dim = n_sites*n_orb

    # create 2 index matrix
    Uij_anti = np.zeros((dim, dim))
    Uij_par = np.zeros((dim, dim))
    Uiijj = np.zeros((dim, dim))
    Uijji = np.zeros((dim, dim))

    for i in range(0, dim):
        for j in range(0, dim):
            # the indices in VASP are switched: U_ijkl ---VASP--> U_ikjl
            Uij_anti[i, j] = uijkl[i, i, j, j]
            Uij_par[i, j] = uijkl[i, i, j, j]-uijkl[i, j, j, i]
            Uiijj[i, j] = uijkl[i, j, i, j]
            Uijji[i, j] = uijkl[i, j, j, i]

    np.set_printoptions(precision=3, suppress=True)

    if out:
        print('reduced U anti-parallel = U_mm\'\^oo\' = U_mm\'mm\' matrix : \n', Uij_anti)
        print('reduced U parallel = U_mm\'\^oo = U_mm\'mm\' - U_mm\'m\'m matrix : \n', Uij_par)
        print('reduced Uijji : \n', Uijji)
        print('reduced Uiijj : \n', Uiijj)

    return Uij_anti, Uiijj, Uijji, Uij_par


def calc_kan_params(uijkl, n_sites, n_orb, out=False):
    '''
    calculates the kanamori interaction parameters from a
    given Uijkl matrix. Follows the procedure given in
    PHYSICAL REVIEW B 86, 165105 (2012) Vaugier,Biermann
    formula 30,31,32

    Parameters
    ----------
    uijkl : numpy array
        4d numpy array of Coulomb tensor
    n_sites: int
        number of different atoms (Wannier centers)
    n_orb : int
        number of orbitals per atom
    out : bool
        verbose mode

    Returns
    -------
    int_params : direct
        kanamori parameters
    '''

    int_params = collections.OrderedDict()

    # calculate intra-orbital U
    U = 0.0
    for i in range(0, n_orb):
        U += uijkl[i, i, i, i]
    U = U/(n_orb)
    int_params['U'] = U

    # calculate the U'
    Uprime = 0.0
    for i in range(0, n_orb):
        for j in range(0, n_orb):
            if i != j:
                Uprime += uijkl[i, i, j, j]
    Uprime = Uprime / (n_orb*(n_orb-1))
    int_params['Uprime'] = Uprime

    # calculate J
    J = 0.0
    for i in range(0, n_orb):
        for j in range(0, n_orb):
            if i != j:
                J += uijkl[i, j, i, j]
    J = J / (n_orb*(n_orb-1))
    int_params['J'] = J

    if out:
        print('U= ', "{:.4f}".format(U))
        print('U\'= ', "{:.4f}".format(Uprime))
        print('J= ', "{:.4f}".format(J))

    return int_params


def fit_kanamori(uijkl, n_orb, switch_jk=False, fit_2=True, fit_3=False, fit_4=True):
    '''
    Fit Kanamori Hamiltonian with scipy to 2,3, and / or 4 parameters

    Parameters
    -----------
    uijkl: np.array (n_orb x n_orb x n_orb x n_orb)
            input four index tensor
    n_orb: int
            number of orbitals
    switch_jk: bool, default=False
            flip two inner indices in input U tensor (for Vasp)
    fit_2: bool, default=True
            fit two parameter form
    fit_3: bool, default=False
            fit three parameter form (U,Up,J=Jc)
    fit_4: bool, default=True
            fit four parameter form

    Returns
    -------
    Uijkl_fit: np.array (n_orb x n_orb x n_orb x n_orb)
            fitted Uijkl tensor
    '''
    from scipy.optimize import minimize

    def minimizer_2params(parameters):
        U, J = parameters

        Uijkl_fit = construct_U_kan(n_orb, U, J)

        return np.sum((uijkl - Uijkl_fit)**2)

    def minimizer_3params(parameters):
        U, J, Up = parameters

        Uijkl_fit = construct_U_kan(n_orb, U, J, Up)

        return np.sum((uijkl - Uijkl_fit)**2)

    def minimizer_4params(parameters):
        U, J, Up, Jc = parameters

        Uijkl_fit = construct_U_kan(n_orb, U, J, Up, Jc)

        return np.sum((uijkl - Uijkl_fit)**2)

    # check if J = JC (Hunds exchange and pair hopping have same amplitude)
    # true for real values wave functions
    if np.max(uijkl.imag) > 0.0:
        print(f"Largest imaginary part of Uijkl: {np.max(uijkl.imag)}. Kanamori Hint assumed to be real valued. Neglecting imag part")
        uijkl = uijkl.real

    if switch_jk:
        uijkl = np.moveaxis(uijkl, 1, 2)

    # fit U, J
    if fit_2:
        initial_guess = (4, 1)
        result = minimize(minimizer_2params, initial_guess)
        U, J = result.x
        Uijkl_fit = construct_U_kan(n_orb, U, J)
        print('Result 2 parameter fit: \nU = {:.4f} eV, J = {:.4f} eV'.format(U, J))
        print(f'optimize error {result.fun:.3e}')
        max_ind = np.unravel_index(np.argmax(np.abs(Uijkl_fit-uijkl), axis=None), Uijkl_fit.shape)
        print(f'U max diff: U{max_ind}= {np.abs(Uijkl_fit-uijkl)[max_ind]:.4e}')

        print('\n-------------------------\n')

    # fit U, J, Up
    if fit_3:
        initial_guess = (4, 1, 2)
        result = minimize(minimizer_3params, initial_guess)
        U, J, Up = result.x
        Uijkl_fit = construct_U_kan(n_orb, U, J, Up)
        print('Result 3 parameter fit: \nU = {:.4f} eV, U\' = {:.4f} eV J = {:.4f} eV'.format(U, Up, J))
        print(f'optimize error {result.fun:.3e}')
        max_ind = np.unravel_index(np.argmax(np.abs(Uijkl_fit-uijkl), axis=None), Uijkl_fit.shape)
        print(f'U max diff: U{max_ind}= {np.abs(Uijkl_fit-uijkl)[max_ind]:.4e}')
        print(f'U=U\'-2J deviation: {U-Up-2*J:.4f}')

        print('\n-------------------------\n')

    if fit_4:
        # fit U, J, Up, Jc
        initial_guess = (4, 1, 2, 1)
        result = minimize(minimizer_4params, initial_guess)
        U, J, Up, Jc = result.x
        Uijkl_fit = construct_U_kan(n_orb, U, Up, J, Jc)
        print('Result 4 parameter fit: \nU = {:.4f} eV, U\' = {:.4f} eV, J = {:.4f} eV, Jc = {:.4f} eV'.format(U, Up, J, Jc))
        print(f'optimize error {result.fun:.3e}')
        print(f'U max diff: U{max_ind}= {np.abs(Uijkl_fit-uijkl)[max_ind]:.4e}')
        max_ind = np.unravel_index(np.argmax(np.abs(Uijkl_fit-uijkl), axis=None), Uijkl_fit.shape)
        print(f'U=U\'-2J deviation: {U-Up-2*J:.4f}')
        print(f'J=Jc deviation: {J-Jc:.4f}')

    return Uijkl_fit


def calc_u_avg_fulld(uijkl, n_sites, n_orb, out=False):
    '''
    calculates the coulomb integrals from a
    given Uijkl matrix for full d shells. Follows the procedure given
    in Pavarini  - 2014 - arXiv - 1411 6906 - julich school U matrix
    page 8 or as done in
    PHYSICAL REVIEW B 86, 165105 (2012) Vaugier,Biermann
    formula 23, 25
    works atm only for full d shell (l=2)

    Returns F0=U, and J=(F2+F4)/14

    Parameters
    ----------
    uijkl : numpy array
        4d numpy array of Coulomb tensor
    n_sites: int
        number of different atoms (Wannier centers)
    n_orb : int
        number of orbitals per atom
    out : bool
        verbose mode

    Returns
    -------
    int_params : direct
        Slater parameters
    '''

    int_params = collections.OrderedDict()
    Uij_anti, Uiijj, Uijji, Uij_par = red_to_2ind(uijkl, n_sites, n_orb, out=out)
    # U_antipar = U_mm'^oo' = U_mm'mm' (Coulomb Int)
    # U_par = U_mm'^oo = U_mm'mm' - U_mm'm'm (for intersite interaction)
    # here we assume cubic harmonics (real harmonics) as basis functions in the order
    # dz2 dxz dyz dx2-y2 dxy
    # triqs basis: basis ordered as (xy,yz,z^2,xz,x^2-y^2)

    # calculate J
    J_cubic = 0.0
    for i in range(0, n_orb):
        for j in range(0, n_orb):
            if i != j:
                J_cubic += Uijji[i, j]
    J_cubic = J_cubic/(20.0)
    # 20 for 2l(2l+1)
    int_params['J_cubic'] = J_cubic

    # conversion from cubic to spherical:
    J = 7.0 * J_cubic / 5.0

    int_params['J'] = J

    # calculate intra-orbital U
    U_0 = 0.0
    for i in range(0, n_orb):
        U_0 += Uij_anti[i, i]
    U_0 = U_0 / float(n_orb)
    int_params['U_0'] = U_0

    # now conversion from cubic to spherical
    U = U_0 - (8.0*J_cubic/5.0)

    int_params['U'] = U

    if out:
        print('cubic U_0= ', "{:.4f}".format(U_0))
        print('cubic J_cubic= ', "{:.4f}".format(J_cubic))
        print('spherical F0=U= ', "{:.4f}".format(U))
        print('spherical J=(F2+f4)/14 = ', "{:.4f}".format(J))

    return int_params


def calculate_interaction_from_averaging(uijkl, n_sites, n_orb, out=False):
    '''
    calculates U,J by averaging directly the Uijkl matrix
    ignoring if tensor is given in spherical or cubic basis.
    The assumption here is that the averaging gives indepentendly
    of the choosen basis (cubic or spherical harmonics) the same results
    if Uijkl is a true Slater matrix.

    Returns F0=U, and J=(F2+F4)/14

    Parameters
    ----------
    uijkl : numpy array
        4d numpy array of Coulomb tensor
    n_sites: int
        number of different atoms (Wannier centers)
    n_orb : int
        number of orbitals per atom
    out : bool
        verbose mode

    Returns
    -------
    U, J: tuple
        Slater parameters
    '''

    l = 2

    Uij_anti, Uiijj, Uijji, Uij_par = red_to_2ind(uijkl, n_sites, n_orb, out=out)

    # Calculates Slater-averaged parameters directly
    U = [None] * n_sites
    J = [None] * n_sites
    for impurity in range(n_sites):
        u_ijij_imp = Uij_anti[impurity*n_orb:(impurity+1)*n_orb, impurity*n_orb:(impurity+1)*n_orb]
        U[impurity] = np.mean(u_ijij_imp)

        u_iijj_imp = Uiijj[impurity*n_orb:(impurity+1)*n_orb, impurity*n_orb:(impurity+1)*n_orb]
        J[impurity] = np.sum(u_iijj_imp) / (2*l*(2*l+1)) - U[impurity] / (2*l)
    U = np.mean(U)
    J = np.mean(J)

    if out:
        print('spherical F0=U= ', "{:.4f}".format(U))
        print('spherical J=(F2+f4)/14 = ', "{:.4f}".format(J))

    return U, J


def fit_slater_fulld(uijkl, n_sites, U_init, J_init, fixed_F4_F2=True):
    '''
    finds best Slater parameters U, J for given Uijkl tensor
    using the triqs U_matrix operator routine
    assumes F4/F2=0.625
    '''

    from triqs.operators.util.U_matrix import U_matrix_slater, reduce_4index_to_2index
    from scipy.optimize import minimize
    # transform U matrix orbital basis ijkl to nmop, note the last two indices need to be switched in the T matrices

    def transformU(U_matrix, T):
        return np.einsum("im,jn,ijkl,lo,kp->mnpo", np.conj(T), np.conj(T), U_matrix, T, T)

    def minimizer(parameters):
        U_int, J_hund = parameters
        Umat_full = U_matrix_slater(l=2, U_int=U_int, J_hund=J_hund, basis='cubic')
        Umat_full = transformU(Umat_full, rot_def_to_w90)

        Umat, Upmat = reduce_4index_to_2index(Umat_full)
        u_iijj_crpa = Uiijj[:5, :5]
        u_iijj_slater = Upmat - Umat
        u_ijij_crpa = Uij_anti[:5, :5]
        u_ijij_slater = Upmat
        return np.sum((u_iijj_crpa - u_iijj_slater)**2 + (u_ijij_crpa - u_ijij_slater)**2)

    def minimizer_radial(parameters):
        F0, F2, F4 = parameters
        Umat_full = U_matrix_slater(l=2, radial_integrals=[F0, F2, F4], basis='cubic')
        Umat_full = transformU(Umat_full, rot_def_to_w90)

        Umat, Upmat = reduce_4index_to_2index(Umat_full)
        u_iijj_crpa = Uiijj[:5, :5]
        u_iijj_slater = Upmat - Umat
        u_ijij_crpa = Uij_anti[:5, :5]
        u_ijij_slater = Upmat
        return np.sum((u_iijj_crpa - u_iijj_slater)**2 + (u_ijij_crpa - u_ijij_slater)**2)

    # rot triqs d basis to w90 default basis!
    # check your order of orbitals assuming:
    # dz2, dxz, dyz, dx2-y2, dxy
    rot_def_to_w90 = np.array([[0, 0, 0, 0, 1],
                               [0, 0, 1, 0, 0],
                               [1, 0, 0, 0, 0],
                               [0, 1, 0, 0, 0],
                               [0, 0, 0, 1, 0]])

    Uij_anti, Uiijj, Uijji, Uij_par = red_to_2ind(uijkl, n_sites, n_orb=5, out=False)

    if fixed_F4_F2:
        result = minimize(minimizer, (U_init, J_init))

        U_int, J_hund = result.x
        print('Final results from fit: U = {:.3f} eV, J = {:.3f} eV'.format(U_int, J_hund))
        print('optimize error', result.fun)
    else:
        # start with 0.63 as guess
        F0 = U_init
        F2 = J_init * 14.0 / (1.0 + 0.63)
        F4 = 0.630 * F2

        initial_guess = (F0, F2, F4)

        print('Initial guess: F0 = {0[0]:.3f} eV, F2 = {0[1]:.3f} eV, F4 = {0[2]:.3f} eV'.format(initial_guess))

        result = minimize(minimizer_radial, initial_guess)
        F0, F2, F4 = result.x
        print('Final results from fit: F0 = {:.3f} eV, F2 = {:.3f} eV, F4 = {:.3f} eV'.format(F0, F2, F4))
        print('(F2+F4)/14 = {:.3f} eV'.format((F2+F4)/14))
        print('F4/F2 = {:.3f} eV'.format(F4/F2))
        print('optimize error', result.fun)
        U_int = F0
        J_hund = (F2+F4)/14

    return U_int, J_hund

# example for a five orbital model
# uijkl=read_uijkl('UIJKL',1,5)
# calc_u_avg_fulld(uijkl,n_sites=1,n_orb=5,out=True)
# print('now via fitting with fixed F4/F2=0.63 ratio')
# fit_slater_fulld(uijkl,1,3,1, fixed_F4_F2 = True)
# print('now directly fitting F0,F2,F4')
# fit_slater_fulld(uijkl,1,3,1, fixed_F4_F2 = False)

# calculate_interaction_from_averaging(uijkl, 1, 5, out=True)


# example for 3 orbital kanamori
# uijkl=read_uijkl('UIJKL',1,3)
# calc_kan_params(uijkl,1,3,out=True)
