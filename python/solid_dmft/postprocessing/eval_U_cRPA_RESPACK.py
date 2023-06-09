#!@TRIQS_PYTHON_EXECUTABLE@

import numpy as np
from itertools import product


class respack_data:
    '''
    respack data class
    '''

    def __init__(self, path, seed):
        self.path = path
        self.seed = seed
        self.freq = None
        self.n_orb = None
        # Uijij
        self.U_R = None
        # Vijij
        self.V_R = None
        # Uijji = Uiijj
        self.J_R = None
        # Vijji = Viijj
        self.X_R = None

        # full Uijkl reconstructed
        self.Uijkl = None
        self.Vijkl = None

        # freq dependent direction Coulomb
        self.Uij_w = None
        self.Jij_w = None
        self.w_mesh = None


def _read_R_file(file):
    '''
    read respack Wmat, Jmat, Vmat, Xmat file format

    Parameters:
    -----------
    file: string
        string to file

    Returns:
    --------
    U_R: dict of np.ndarray
        keys are tuples of 3d integers representing the R vectors
        values are n_orb x n_orb np.ndarray complex
    n_orb: int
        number of orbitals
    '''
    # read X_R from file
    with open(file, 'r') as fd:
        # eliminate header
        fd.readline()
        fd.readline()
        fd.readline()

        # get number of R vectors
        r_vec_max = [abs(int(i)) for i in fd.readline().strip().split()]
        assert len(r_vec_max) == 3

        # determine number of orbitals
        while True:
            line = fd.readline().strip().split()
            if not line:
                break
            else:
                n_orb = int(line[0])

    # open again and read whole file
    U_R = {}
    with open(file, 'r') as fd:
        # eliminate header
        fd.readline()
        fd.readline()
        fd.readline()

        for x, y, z in product(range(-r_vec_max[0], r_vec_max[0]+1),
                               range(-r_vec_max[1], r_vec_max[1]+1),
                               range(-r_vec_max[2], r_vec_max[2]+1)
                               ):
            fd.readline()  # remove rvec line
            U_R[tuple([x, y, z])] = np.zeros((n_orb, n_orb), dtype=complex)
            for i, j in product(range(n_orb), range(n_orb)):
                line = fd.readline().strip().split()
                U_R[tuple([x, y, z])][i, j] = float(line[2])+1j*float(line[3])
            fd.readline()  # remove empty line before next block

    return U_R, n_orb


def _read_freq_int(path, n_orb, U_or_J, w_or_iw='w'):
    '''
    read frequency dependent files from disk

    Parameters:
    -----------
    path: string
        path to respack calculations
    n_orb: int
        number of orbitals
    U_or_J: string
        pass either U or J for reading U or reading J
    w_or_iw: string, optional default='w'
        read either real frequency axis or Matsubara axis results

    Returns:
    --------
    Uij_w: np.ndarray, dtype=complex shape=(n_w, n_orb, n_orb)
        direct freq dependent Coulomb integrals between orbitals
    w_mesh: np.ndarry, shape=(n_w)
        frequency mesh of Uij_w tensor
    '''

    assert U_or_J == 'U' or U_or_J == 'J'

    if U_or_J == 'U':
        file = path+'/dir-intW/dat.UvsE.'
    else:
        file = path+'/dir-intJ/dat.JvsE.'

    # read first w_mesh
    if w_or_iw == 'w':
        w_mesh = np.loadtxt(file+'001-001')[:, 0]
    else:
        w_mesh = np.loadtxt(file+'001-001')[:, 1]

    Uij_w = np.zeros((w_mesh.shape[0], n_orb, n_orb), dtype=complex)

    for i in range(0, n_orb):
        for j in range(i, n_orb):  # only the upper triangle part is calculated by RESPACK
            temp_u = np.loadtxt(file+str(i+1).zfill(3)+'-'+str(j+1).zfill(3))
            Uij_w[:, i, j] = temp_u[:, 2] + 1j*temp_u[:, 3]
            if not i == j:  # set the lower triangle with the complex conj
                Uij_w[:, j, i] = temp_u[:, 2] + - 1j*temp_u[:, 3]

    return Uij_w, w_mesh


def read_interaction(seed, path='./'):
    '''
    Parameters:
    -----------
    seed: string
        seed of QE / w90 file names
    path: string
        path to respack calculations

    Returns:
    --------
    res: respack data class
    '''

    res = respack_data(seed, path)

    # read w3d.out file for frequency info
    with open(path+'/'+seed+'.w3d.out', 'r') as w3d:
        lines = w3d.readlines()
        for line in lines:
            if 'CALC_IFREQ=' in line:
                res.freq = int(line.replace(' ', '').strip().split('=')[1].split(',')[0])
                break

    # read R dependent matrices and IFREQ
    res.U_R, res.n_orb = _read_R_file(file=path+'/dir-intW/dat.Wmat')
    res.V_R, _ = _read_R_file(file=path+'/dir-intW/dat.Vmat')
    res.J_R, _ = _read_R_file(file=path+'/dir-intJ/dat.Jmat')
    res.X_R, _ = _read_R_file(file=path+'/dir-intJ/dat.Xmat')

    # create R=0 matrices for user convenience
    res.Uijij = res.U_R[(0, 0, 0)]
    res.Uijji = res.J_R[(0, 0, 0)]
    res.Vijij = res.V_R[(0, 0, 0)]
    res.Vijji = res.X_R[(0, 0, 0)]

    # reconstruct full Uijkl tensor assuming Uijji == Uiijj
    res.Uijkl = construct_Uijkl(res.U_R[(0, 0, 0)], res.J_R[(0, 0, 0)])
    res.Vijkl = construct_Uijkl(res.V_R[(0, 0, 0)], res.X_R[(0, 0, 0)])

    # read freq dependent results
    res.Uij_w, res.w_mesh = _read_freq_int(path, res.n_orb, U_or_J='U')
    res.Jij_w, res.w_mesh = _read_freq_int(path, res.n_orb, U_or_J='J')

    return res


def construct_Uijkl(Uijij, Uiijj):
    '''
    construct full 4 index Uijkl tensor from respack data
    assuming Uijji = Uiijj

    ----------
    Uijij: np.ndarray
        Uijij matrix
    Uiijj: np.ndarray
        Uiijj matrix

    -------
    uijkl : numpy array
        uijkl Coulomb tensor

    '''

    n_orb = Uijij.shape[0]
    orb_range = range(0, n_orb)
    Uijkl = np.zeros((n_orb, n_orb, n_orb, n_orb), dtype=complex)

    for i, j, k, l in product(orb_range, orb_range, orb_range, orb_range):
        if i == j == k == l:  # Uiiii
            Uijkl[i, j, k, l] = Uijij[i, j]
        elif i == k and j == l:  # Uijij
            Uijkl[i, j, k, l] = Uijij[i, j]
        elif i == l and j == k:  # Uijji
            Uijkl[i, j, k, l] = Uiijj[i, j]
        elif i == j and k == l:  # Uiijj
            Uijkl[i, j, k, l] = Uiijj[i, k]
    return Uijkl


def fit_slater_fulld(u_ijij_crpa, u_ijji_crpa, U_init, J_init, fixed_F4_F2=True):
    '''
    finds best Slater parameters U, J for given Uijij and Uijji matrices
    using the triqs U_matrix operator routine assumes F4/F2=0.625

    Parameters:
    -----------
    u_ijij_crpa: np.ndarray of shape (5, 5)
        Uijij matrix
    u_ijji_crpa: np.ndarray of shape (5, 5)
        Uijji matrix
    U_init: float
        inital value of U for optimization
    J_init: float
        inital value of J for optimization
    fixed_F4_F2: bool, optional default=True
        fix F4/F2 ratio to 0.625
    Returns:
    --------
    U_int: float
        averaged U value
    J_hund: float
        averaged J value
    '''

    from triqs.operators.util.U_matrix import U_matrix_slater, spherical_to_cubic, reduce_4index_to_2index, transform_U_matrix
    from scipy.optimize import minimize

    # input checks
    assert u_ijij_crpa.shape == (5, 5), 'fit slater only implemented for full d shell (5 orbitals)'
    assert u_ijji_crpa.shape == (5, 5), 'fit slater only implemented for full d shell (5 orbitals)'

    def minimizer(parameters):
        U_int, J_hund = parameters
        Umat_full = U_matrix_slater(l=2, U_int=U_int, J_hund=J_hund, basis='spherical')
        Umat_full = transform_U_matrix(Umat_full, T)

        Umat, u_ijij_slater = reduce_4index_to_2index(Umat_full)
        u_iijj_slater = u_ijij_slater - Umat
        return np.sum((u_ijji_crpa - u_iijj_slater)**2 + (u_ijij_crpa - u_ijij_slater)**2)

    def minimizer_radial(parameters):
        F0, F2, F4 = parameters
        Umat_full = U_matrix_slater(l=2, radial_integrals=[F0, F2, F4], basis='spherical')
        Umat_full = transform_U_matrix(Umat_full, T)

        Umat, u_ijij_slater = reduce_4index_to_2index(Umat_full)
        u_iijj_slater = u_ijij_slater - Umat
        return np.sum((u_ijji_crpa - u_iijj_slater)**2 + (u_ijij_crpa - u_ijij_slater)**2)

    # transformation matrix from spherical to w90 basis
    T = spherical_to_cubic(l=2, convention='wannier90')

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
