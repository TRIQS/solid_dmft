from typing import NamedTuple, Any

import numpy as np
from pydlr import kernel, dlr  # TODO: switch from pydlr backend to TRIQS dlr backend
from scipy.optimize import minimize, LinearConstraint, NonlinearConstraint

from triqs.gf import BlockGf, Gf, GfImTime, GfImFreq, MeshImFreq, MeshImTime

is_block_gf = lambda x : isinstance(x, BlockGf)
is_array = lambda x : isinstance(x, np.ndarray)

class CallbackResult(NamedTuple):
    x  : Any      = None
    sigma : Any   = None
    residual : Any = None

class MinimizerResult(NamedTuple):
    scipy_sol : Any = None
    sig_xaa : Any  = None
    g_xaa   : Any  = None
    g0_xaa  : Any  = None
    callback : CallbackResult = None

class SolverResult(NamedTuple):
    Sigma_iw  : Any    = None
    G_tau     : Any    = None
    G0_tau    : Any    = None
    Sigma_moments : Any = None
    minimizer : MinimizerResult     = None

class Callback:

    def __init__(self, eval_func, diff_func, method):
        self._history   = [] # stores (x⃗, Σiνₖ, G-G₀-G₀ΣG) in a CallbackResult
        self._eval_func = eval_func
        self._diff_func = diff_func
        self._method     = method

    def __call__(self, *args):
        x = args[0]
        self._history.append(CallbackResult(x=x, sigma=self._eval_func(x),
                                             residual=self._diff_func(x)))
    @property
    def history(self): return self._history


class Symmetrizer:

    def __init__(self, nx, no):
        self.N = (no*(no-1))//2
        self.nx, self.no = nx, no
        self.diag_idxs = np.arange(self.no)
        self.triu_idxs = np.triu_indices(no, k=1)
        self.tril_idxs = np.tril_indices(no, k=-1)
    
    def get_x_d(self, g_xaa):
        x_d = g_xaa[:, self.diag_idxs, self.diag_idxs].flatten()
        return x_d

    def set_x_d(self, g_xaa, x_d):
        g_xaa[:, self.diag_idxs, self.diag_idxs] = x_d.reshape((self.nx, self.no))
        return g_xaa

    def get_x_u(self, g_xaa):
        x_u = g_xaa[:, self.triu_idxs[0], self.triu_idxs[1]].flatten()
        return x_u

    def set_x_u(self, g_xaa, x_u):
        g_xaa[:, self.triu_idxs[0], self.triu_idxs[1]] = x_u.reshape((self.nx, self.N))
        g_xaa[:, self.tril_idxs[0], self.tril_idxs[1]] = g_xaa[:, self.triu_idxs[0], self.triu_idxs[1]].conj()
        return g_xaa
    
    def get_x_l(self, g_xaa):
        x_l = g_xaa[:, self.tril_idxs[0], self.tril_idxs[1]].flatten()
        return x_l
    
    def set_x_l(self, g_xaa, x_l):
        g_xaa[:, self.tril_idxs[0], self.tril_idxs[1]] = x_l.reshape((self.nx, self.N))
        return g_xaa
        
    def get_diag_indices(self): return self.diag_idxs
    def get_triu_indices(self): return self.triu_idxs



class Dyson:

    def __init__(self, lamb, 
                       eps=1e-9,                
                       method='trust-constr', 
                       options=dict(maxiter=10000),
                       verbose=True,
                       **kwargs
                       ):

        self.lamb    = lamb                                        # dlr lambda
        self.eps     = eps                                         # dlr epsilon
        self._dlr    = dlr(lamb=self.lamb, eps=self.eps, **kwargs) # dlr class instance
        self.Mkl     = self._compute_mkl()                         # compute residual matrix
        self.method  = method                                      # scipy minimize method
        self.options = options                                     # minimize options
        self.verbose = verbose

        # if verbose add display to minimize options
        if self.verbose: self.options['disp'] = True
        else: self.options['disp'] = False

    def __len__(self): return len(self._dlr)

    def __repr__(self):
        out = '-'*20 + ' Dyson solver ' + '-'*20
        out += '\nΛ = {} (Λ ∼ ωmax*β))'.format(self.lamb)
        out += '\nε = {:1.0e}'.format(self.eps)
        out += '\ndlr basis size = {}'.format(len(self))
        out += '\nscipy info: '
        out += '\n\tmethod = {}'.format(self.method)
        for key, val in self.options.items():
            out += '\n\t{}     =    {}'.format(key, val)
        out += '\n'+'-'*20 + '--------------' + '-'*20
        return out

    __str__ = __repr__

    # wrapping functions to the DLR class
    # abstract away the underlying DLR class functions
    @property
    def rank(self): return self._dlr.rank

    @property
    def dlr_nodes(self): return self._dlr.dlrrf

    def get_tau(self, beta): return self._dlr.get_tau(beta)
    def get_iom(self, beta): return self._dlr.get_matsubara_frequencies(beta)
    
    # dlr <-> Matsubara
    def dlr_from_iom(self, g_iwaa, beta): return self._dlr.dlr_from_matsubara(g_iwaa, beta)
    def iom_from_dlr(self, g_xaa,beta): return self._dlr.matsubara_from_dlr(g_xaa, beta)
    def fit_dlr_from_iom(self, mesh, g_iwaa, beta): return self._dlr.lstsq_dlr_from_matsubara(mesh, g_iwaa, beta)
    def eval_dlr_iom(self, g_xaa, mesh, beta): return self._dlr.eval_dlr_freq(g_xaa, mesh, beta)

    # dlr <-> tau
    def dlr_from_tau(self, g_xaa): return self._dlr.dlr_from_tau(g_xaa)
    def fit_dlr_from_tau(self, mesh, g_iaa, beta): return self._dlr.lstsq_dlr_from_tau(mesh, g_iaa, beta)
    def eval_dlr_tau(self, g_xaa, mesh, beta): return self._dlr.eval_dlr_tau(g_xaa, mesh, beta)


    # precomputes Mkl for L2 norm in τ
    def _compute_mkl(self):
        Mkl = np.zeros((len(self), len(self)), dtype=np.float128)
        for iwk, wk in enumerate(self.dlr_nodes):
            for iwl, wl in enumerate(self.dlr_nodes):
                K0wk, Kbwk = kernel(np.array([0.,1.]), np.array([wk]))
                K0wl, Kbwl = kernel(np.array([0.,1.]), np.array([wl]))
                if np.fabs(wk+wl) < 1e-13: 
                    Mkl[iwk,iwl] = K0wk*K0wl
                else: 
                    Mkl[iwk, iwl] = (K0wk*K0wl - Kbwk*Kbwl)
                    Mkl[iwk, iwl] /= (wk+wl)
        return Mkl


    # run the scipy minimization
    def _minimize_dyson_equation(self,
                                 g_iaa,         # G data
                                 g0_iaa,        # G0 data
                                 beta,          # inverse temperature
                                 sigma_moments, # high-freq moments of Σ
                                 tau=None,      # tau mesh
                                ):

        
        # fold and unfold complex numbers
        def merge_re_im(x):
            x_d, x_u = x[:nx*no], x[nx*no:]
            re, im = np.split(x_u, 2)
            x_u = re + 1.j * im
            return x_d, x_u

        def split_re_im(x_d, x_u):
            return np.concatenate((
                np.array(x_d.real, dtype=float),
                np.array(x_u.real, dtype=float),
                np.array(x_u.imag, dtype=float)))
                                       
        # convert dlr coefficients vector to back to complex sigma
        def sig_from_x(x):
            x_d, x_off = x[:2*nx*no], x[2*nx*no:]
            x_u, x_l = np.split(x_off, 2)
            sig = np.zeros((nx, no, no), dtype=dtype)
            for func, x in zip((sym.set_x_d, sym.set_x_u, sym.set_x_l), (x_d, x_u, x_l)):
                re, im = np.split(x,2)
                x = re + 1.j*im
                func(sig, x)
            return sig

        # convert complex sigma to dlr coefficients vector 
        def x_from_sig(sig):
            x_d, x_u, x_l = sym.get_x_d(sig), sym.get_x_u(sig), sym.get_x_l(sig);
            x = np.concatenate((np.array(x_d.real, dtype=float),
                            np.array(x_d.imag, dtype=float),
                            np.array(x_u.real, dtype=float),
                            np.array(x_u.imag, dtype=float),
                            np.array(x_l.real, dtype=float),
                            np.array(x_l.imag, dtype=float)
                          ))
        
            return x

        # matrix to vector transform
        def mat_vec(mat):
            v_d, v_u, v_l = sym.get_x_d(mat[None,...]), sym.get_x_u(mat[None,...]), sym.get_x_l(mat[None,...])
            return np.concatenate((
                np.array(v_d.real, dtype=float),
                np.array(v_d.imag, dtype=float),
                np.array(v_u.real, dtype=float),
                np.array(v_u.imag, dtype=float),
                np.array(v_l.real, dtype=float),
                np.array(v_l.imag, dtype=float)
            ))
            
        # constraint condition: -∑σk =  Σ_1
        def constraint_func(x):
            sig = self.dlr_from_iom(sig_from_x(x), beta)
            mat = -sig.sum(axis=0) 
            vec = mat_vec(mat)
            return vec
        

        # target function  
        def dyson_difference(x):
            sig = sig_from_x(x)
            sig_iwaa = sig + sig_infty
            #  G - G0 - G0*Σ*G = 0 done on the DLR nodes
            r_iwaa = g_iwaa - g0_iwaa - g0_iwaa@sig_iwaa@g_iwaa
            # compute DLR of rk_iwaa
            r_xaa = self.fit_dlr_from_iom(freq, r_iwaa, beta)
            # ||R||^2 = r^T @ M @ r
            R2 = np.einsum('mnk, kl, lnm->nm', r_xaa.T.conj(), self.Mkl, r_xaa).flatten()
            # the Frobeinus norm
            return np.sqrt(np.sum(R2)).real

        assert g_iaa.shape[1:] == g0_iaa.shape[1:] == sigma_moments.shape[1:], "number of orbs inconsistent across G, G0, and moments"
        
        nx = len(self)
        ni, no, _ = g_iaa.shape
        shape_xaa = (nx, no, no)
        N = (no*(no-1))//2
        dtype = complex

        # self energy <-> vector conversion
        sym = Symmetrizer(nx, no)

        # constraint
        sig_infty, sigma_1 = sigma_moments[0], sigma_moments[1]
        bound = mat_vec(sigma_1)
        constraints = (NonlinearConstraint(constraint_func, bound, bound))
        
        # frequency grid points for evaluating dyson equation
        freq = self.get_iom(beta)
        
        # dlr fit to G and G0 
        if tau is None: g_xaa, g0_xaa = self.dlr_from_tau(g_iaa), self.dlr_from_tau(g0_iaa)
        else: g_xaa, g0_xaa  = self.fit_dlr_from_tau(tau, g_iaa, beta), self.fit_dlr_from_tau(tau, g0_iaa, beta)
        
        if self.verbose:
            eval_tau = tau if tau is not None else self.get_tau(beta)
            g, g0 = self.eval_dlr_tau(g_xaa, eval_tau, beta), self.eval_dlr_tau(g0_xaa, eval_tau,  beta)
            print('initial DLR fits to G(τ) and G0(τ)')
            print(f'max|G(τ) - Gdlr(τ)| = {np.max(np.abs(g-g_iaa)):.6e}')
            print(f'max|G0(τ) - G0dlr(τ)| = {np.max(np.abs(g0-g0_iaa)):.6e}')
        
        # compute and obtain initial Σ
        g_iwaa, g0_iwaa = self.iom_from_dlr(g_xaa, beta), self.iom_from_dlr(g0_xaa, beta)
        
        # the DLR representable part of the self-energy
        # initial guess for DLR Sigma only use of Dyson equation!
        sig0_iwaa = np.linalg.inv(g0_iwaa)-np.linalg.inv(g_iwaa)-sig_infty
       
        # optimize Σ(iν)
        x_init = x_from_sig(sig0_iwaa)
        
        cb = Callback(eval_func= lambda x : sig_from_x(x)+sig_infty, diff_func=dyson_difference, method=self.method)

        solution = minimize(dyson_difference, 
                            x_init,
                            method=self.method,
                            constraints=constraints,
                            options=self.options,
                            callback = cb
                           )
        
        if self.verbose: print(solution.message)
        if not solution.success: print('[WARNING] Minimization did not converge! Please proceed with caution!')
        
        sig_iwaa = sig_from_x(solution.x)
        sig_xaa = self.dlr_from_iom(sig_iwaa, beta)
            
        if self.verbose: print(f'Σ1 constraint diff: {np.max(np.abs(-sig_xaa.sum(axis=0)-sigma_1)):.4e}')

        result = MinimizerResult(solution,
                                  sig_xaa,
                                  g_xaa,
                                  g0_xaa,
                                  cb.history
                                 )
        return result

    # main solve function
    def solve(self, Sigma_iw=None, 
                    G_tau=None, 
                    G0_tau=None, 
                    Sigma_moments=None, 
                    beta=None, 
                    om_mesh=None, 
                    tau_mesh=None
                    ):

        result = None
        
        # we are working with a TRIQS Green's function/ Block Green's function object
        if all(list(map(is_block_gf, [G_tau, G0_tau]))):

            beta = G_tau.mesh.beta
            tau  = np.array([float(x) for x in G_tau.mesh]) if tau_mesh is None else tau_mesh

            #TODO: Sigma_iw shouldn't be required for BlockGf option
            if Sigma_iw is not None:
                Sigma_iw_fit = Sigma_iw.copy()
                om_mesh = np.array([complex(x) for x in Sigma_iw_fit.mesh])

            dlr_results = {}
            for block, sig in Sigma_iw_fit:

                dlr_results[block] =  self._minimize_dyson_equation(G_tau[block].data,
                                                                    G0_tau[block].data,
                                                                    beta,
                                                                    Sigma_moments[block],
                                                                    tau=tau
                                                              )
                # DLR representable part of self energy
                Sigma_iw_fit[block].data[:] = self.eval_dlr_iom(dlr_results[block].sig_xaa, om_mesh, beta)
                # add constant back to obtain true self energy
                Sigma_iw_fit[block].data[:] +=  Sigma_moments[block][0]

        # our Green's functions are just numpy arrays
        elif all(list(map(is_array, [G_tau, G0_tau]))):

            assert beta is not None, "must provide a beta!"

            dlr_results = self._minimize_dyson_equation(G_tau,
                                                        G0_tau,
                                                        beta,
                                                        Sigma_moments,
                                                        tau=tau_mesh
                                                        ) 

            if om_mesh is None:
                Sigma_iw_fit = self.iom_from_dlr(dlr_results.sig_xaa, beta)
                Sigma_iw_fit += Sigma_moments[0]
            else:
                Sigma_iw_fit = self.eval_dlr_iom(dlr_results.sig_xaa,om_mesh,beta) 
                Sigma_iw_fit += Sigma_moments[0]

        else: raise ValueError

        result = SolverResult(Sigma_iw_fit,
                               G_tau,
                               G0_tau,
                               Sigma_moments,
                               dlr_results
                        )

        return result
