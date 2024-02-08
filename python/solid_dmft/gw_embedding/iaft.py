import numpy as np
import sparse_ir

"""
Fourier transform on the imaginary axis based on IR basis and the sparse sampling technique.  
"""


class IAFT(object):
    """
    Driver for FT on the imaginary axis.
    Given inverse temperature, lambda and precision, the IAFT class evaluate the corresponding
    IR basis and sparse sampling points on-the-fly.

    Dependency:
        sparse-ir with xprec supports (https://sparse-ir.readthedocs.io/en/latest/)
        To install sparse-ir with xprec supports: "pip install sparse-ir[xprex]".

    Attributes:
    beta: float
        Inverse temperature (a.u.)
    lmbda: float
        Dimensionless lambda parameter for constructing the IR basis
    prec: float
        Precision for IR basis
    bases: sparse-ir.FiniteTempBasisSet
        IR basis instance
    tau_mesh_f: numpy.ndarray(dim=1)
        Fermionic tau sampling points
    tau_mesh_b: numpy.ndarray(dim=1)
        Bosonic tau sampling points
    wn_mesh_f: numpy.ndarray(dim=1)
        Fermionic Matsubara "indices" sampling points. NOT PHYSICAL FREQUENCIES.
        Physical Matsubara frequencies are wn_mesh_f * numpy.pi / beta
    wn_mesh_b: numpy.ndarray(dim=1)
        Bosonic Matsubara "indices" sampling points. NOT PHYSICAL FREQUENCIES.
        Physical Matsubara frequencies are wn_mesh_f * numpy.pi / beta
    nt_f: int
        Number of fermionic tau sampling points
    nt_b: int
        Number of bosonic tau sampling points
    nw_f: int
        Number of fermionic frequency sampling points
    nw_b: int
        Number of bosonic frequency sampling points
    """
    def __init__(self, beta: float, lmbda: float, prec: float = 1e-15):
        """
        :param beta: float
            Inverse temperature (a.u.)
        :param lmbda: float
            Lambda parameter for constructing IR basis.
        :param prec: float
            Precision for IR basis
        """
        self.beta = beta
        self.lmbda = lmbda
        self.prec = prec
        self.wmax = lmbda / beta
        self.statisics = {'f', 'b'}

        self.bases = sparse_ir.FiniteTempBasisSet(beta=beta, wmax=self.wmax, eps=prec)
        self.tau_mesh_f = self.bases.smpl_tau_f.sampling_points
        self.tau_mesh_b = self.bases.smpl_tau_b.sampling_points
        self._wn_mesh_f = self.bases.smpl_wn_f.sampling_points
        self._wn_mesh_b = self.bases.smpl_wn_b.sampling_points
        self.nt_f, self.nw_f = self.tau_mesh_f.shape[0], self._wn_mesh_f.shape[0]
        self.nt_b, self.nw_b = self.tau_mesh_b.shape[0], self._wn_mesh_b.shape[0]

        Ttl_ff = self.bases.basis_f.u(self.tau_mesh_f).T
        Twl_ff = self.bases.basis_f.uhat(self._wn_mesh_f).T
        Ttl_bb = self.bases.basis_b.u(self.tau_mesh_b).T
        Twl_bb = self.bases.basis_b.uhat(self._wn_mesh_b).T

        self.Tlt_ff = np.linalg.pinv(Ttl_ff)
        self.Tlt_bb = np.linalg.pinv(Ttl_bb)
        self.Tlw_ff = np.linalg.pinv(Twl_ff)
        self.Tlw_bb = np.linalg.pinv(Twl_bb)

        # Ttw_ff = Ttl_ff * [Twl_ff]^{-1}
        self.Ttw_ff = np.dot(Ttl_ff, self.Tlw_ff)
        self.Twt_ff = np.dot(Twl_ff, self.Tlt_ff)
        self.Ttw_bb = np.dot(Ttl_bb, self.Tlw_bb)
        self.Twt_bb = np.dot(Twl_bb, self.Tlt_bb)

        print(self)

    def __str__(self):
        return "*******************************\n" \
               "Imaginary-Axis Fourier Transform based on IR basis and sparse-sampling:\n" \
               "*******************************\n" \
               "    - precision = {}\n" \
               "    - beta = {}\n" \
               "    - lambda = {}\n" \
               "    - nt_f, nw_f = {}, {}\n" \
               "    - nt_b, nw_b = {}, {}\n" \
               "*******************************".format(self.prec, self.beta, self.lmbda, self.nt_f, self.nw_f,
                                                        self.nt_b, self.nw_b)

    def wn_mesh(self, stats: str, ir_notation: bool = True):
        """
        Return Matsubara frequency indices.
        :param stats: str
            statistics: 'f' for fermions and 'b' for bosons
        :param ir_notation: bool
            Whether wn_mesh_interp is in sparse_ir notation where iwn = n*pi/beta for both fermions and bosons.
            Otherwise, iwn = (2n+1)*pi/beta  for fermions and 2n*pi/beta for bosons.

        :return: numpy.ndarray(dim=1)
            Matsubara frequency indices
        """
        if stats not in self.statisics:
            raise ValueError("Unknown statistics '{}'. "
                             "Acceptable options are 'f' for fermion and 'b' for bosons.".format(stats))
        wn_mesh = np.array(self._wn_mesh_f, dtype=int) if stats == 'f' else np.array(self._wn_mesh_b, dtype=int)
        if not ir_notation:
            wn_mesh = (wn_mesh-1)//2 if stats == 'f' else wn_mesh//2
        return wn_mesh

    def tau_to_w(self, Ot, stats: str):
        """
        Fourier transform from imaginary-time axis to Matsubara-frequency axis
        :param Ot: numpy.ndarray
            imaginary-time object with dimensions (nts, ...)
        :param stats: str
            statistics: 'f' for fermions and 'b' for bosons

        :return: numpy.ndarray
            Matsubara-frequency object with dimensions (nw, ...)
        """
        if stats not in self.statisics:
            raise ValueError("Unknown statistics '{}'. "
                             "Acceptable options are 'f' for fermion and 'b' for bosons.".format(stats))
        Twt = self.Twt_ff if stats == 'f' else self.Twt_bb
        if Ot.shape[0] != Twt.shape[1]:
            raise ValueError(
                "tau_to_w: Number of tau points are inconsistent: {} and {}".format(Ot.shape[0], Twt.shape[1]))

        Ot_shape = Ot.shape
        Ot = Ot.reshape(Ot.shape[0], -1)
        Ow = np.dot(Twt, Ot)

        Ot = Ot.reshape(Ot_shape)
        Ow = Ow.reshape((Twt.shape[0],) + Ot_shape[1:])
        return Ow

    def w_to_tau(self, Ow, stats):
        """
        Fourier transform from Matsubara-frequency axis to imaginary-time axis.

        :param Ow: numpy.ndarray
            Matsubara-frequency object with dimensions (nw, ...)
        :param stats: str
            statistics, 'f' for fermions and 'b' for bosons

        :return: numpy.ndarray
            Imaginary-time object with dimensions (nt, ...)
        """
        if stats not in self.statisics:
            raise ValueError("Unknown statistics '{}'. "
                             "Acceptable options are 'f' for fermion and 'b' for bosons.".format(stats))
        Ttw = self.Ttw_ff if stats == 'f' else self.Ttw_bb
        if Ow.shape[0] != Ttw.shape[1]:
            raise ValueError(
                "w_to_tau: Number of w points are inconsistent: {} and {}".format(Ow.shape[0], Ttw.shape[1]))

        Ow_shape = Ow.shape
        Ow = Ow.reshape(Ow.shape[0], -1)
        Ot = np.dot(Ttw, Ow)

        Ow = Ow.reshape(Ow_shape)
        Ot = Ot.reshape((Ttw.shape[0],) + Ow_shape[1:])
        return Ot

    def w_interpolate(self, Ow, wn_mesh_interp, stats: str, ir_notation: bool = True):
        """
        Interpolate a dynamic object to arbitrary points on the Matsubara axis.

        :param Ow: numpy.ndarray
            Dynamic object on the Matsubara sampling points, self.wn_mesh.
        :param wn_mesh_interp: numpy.ndarray(dim=1, dtype=int)
            Target frequencies "INDICES".
            The physical Matsubara frequencies are wn_mesh_interp * pi/beta.
        :param stats: str
            Statistics, 'f' for fermions and 'b' for bosons.
        :param ir_notation: bool
            Whether wn_mesh_interp is in sparse_ir notation where iwn = n*pi/beta for both fermions and bosons.
            Otherwise, iwn = (2n+1)*pi/beta  for fermions and 2n*pi/beta for bosons.

        :return: numpy.ndarray
            Matsubara-frequency object with dimensions (nw_interp, ...)
        """
        if stats not in self.statisics:
            raise ValueError("Unknown statistics '{}'. "
                             "Acceptable options are 'f' for fermion and 'b' for bosons.".format(stats))
        if ir_notation:
            wn_indices = wn_mesh_interp
        else:
            wn_indices = np.array([2*n+1 if stats == 'f' else 2*n for n in wn_mesh_interp], dtype=int)
        Tlw = self.Tlw_ff if stats == 'f' else self.Tlw_bb
        if Ow.shape[0] != Tlw.shape[1]:
            raise ValueError(
                "w_interpolate: Number of w points are inconsistent: {} and {}".format(Ow.shape[0], Tlw.shape[1]))

        Twl_interp = self.bases.basis_f.uhat(wn_indices).T if stats == 'f' else self.bases.basis_b.uhat(wn_indices).T
        Tww = np.dot(Twl_interp, Tlw)

        Ow_shape = Ow.shape
        Ow = Ow.reshape(Ow.shape[0], -1)
        Ow_interp = np.dot(Tww, Ow)

        Ow = Ow.reshape(Ow_shape)
        Ow_interp = Ow_interp.reshape((wn_indices.shape[0],) + Ow_shape[1:])
        return Ow_interp

    def tau_interpolate(self, Ot, tau_mesh_interp, stats: str):
        """
         Interpolate a dynamic object to arbitrary points on the imaginary-time axis.

        :param Ot: numpy.ndarray
            Dynamic object on the imaginary-time sampling points, self.tau_mesh.
        :param tau_mesh_interp: numpy.ndarray(dim=1, dtype=float)
            Target tau points.
        :param stats: str
            Statistics, 'f' for fermions and 'b' for bosons

        :return: numpy.ndarray
            Imaginary-time object with dimensions (nt_interp, ...)
        """
        if stats not in self.statisics:
            raise ValueError("Unknown statistics '{}'. "
                             "Acceptable options are 'f' for fermion and 'b' for bosons.".format(stats))
        Tlt = self.Tlt_ff if stats == 'f' else self.Tlt_bb
        if Ot.shape[0] != Tlt.shape[1]:
            raise ValueError(
                "t_interpolate: Number of tau points are inconsistent: {} and {}".format(Ot.shape[0], Tlt.shape[1]))

        Ttl_interp = self.bases.basis_f.u(tau_mesh_interp).T if stats == 'f' else self.bases.basis_b.u(tau_mesh_interp).T
        Ttt = np.dot(Ttl_interp, Tlt)

        Ot_shape = Ot.shape
        Ot = Ot.reshape(Ot.shape[0], -1)
        Ot_interp = np.dot(Ttt, Ot)

        Ot = Ot.reshape(Ot_shape)
        Ot_interp = Ot_interp.reshape((tau_mesh_interp.shape[0],) + Ot_shape[1:])
        return Ot_interp


if __name__ == '__main__':
    # Initialize IAFT object for given inverse temperature, lambda and precision
    ft = IAFT(1000, 1e4, 1e-6)

    print(ft.wn_mesh('f', True))

    Gt = np.zeros((ft.nt_f, 2, 2, 2))
    Gw = ft.tau_to_w(Gt, 'f')
    print(Gw.shape)

    # Interpolate to arbitrary tau point
    tau_interp = np.array([0.0, ft.beta])
    Gt_interp = ft.tau_interpolate(Gt, tau_interp, 'f')
    print(Gt_interp.shape)

    # wn in spare_ir notation
    w_interp = np.array([-1,1,3,5], dtype=int)
    Gw_interp = ft.w_interpolate(Gw, w_interp, 'f', True)
    print(Gw_interp.shape)

    # wn in physical notation
    w_interp = np.array([-1,0,1,2,3,4], dtype=int)
    Gw_interp = ft.w_interpolate(Gw, w_interp, 'f', False)
    print(Gw_interp.shape)

    Gt2 = ft.w_to_tau(Gw, 'f')
    print(Gt2.shape)
