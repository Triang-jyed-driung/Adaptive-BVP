import numpy as np
from math import *
from numpy.linalg import lstsq, pinv
import numpy.polynomial.chebyshev as cheb
import random
from matplotlib import pyplot as plt

inv = pinv

def solve(A, b, residuals=False, rank=False, svd=False):
    x, _residuals, _rank, _svd = lstsq(A, b, rcond=None)
    ret = [x]
    if residuals:
        ret += [_residuals]
    if rank:
        ret += [_rank]
    if svd:
        ret += [_svd]
    return (ret[0] if len(ret) == 1 else ret)

K = 16
tau = cheb.chebpts1(K)
# T: chebyshev interpolation matrix (K by K)
# transforms coefficients to points
T = cheb.chebvander(tau, K-1)
TK1 = cheb.chebvander(tau, K)
# T_inv: get coefficients from points
T_inv = inv(T)
# I: Chebyshev integration matrix (K+1 by K)
Il = cheb.chebint(np.eye(K), m=1, k=0, lbnd=-1, axis=0)
IlK = Il[0:K]
Ir = -cheb.chebint(np.eye(K), m=1, k=0, lbnd=1, axis=0)
IrK = Ir[0:K]
cheb_integ = np.array([(2/(1-m*m) if m % 2 == 0 else 0) for m in range(K)])
cheb_quad = cheb_integ @ T_inv

# (K, K+1) @ (K+1, K) @ (K, K)
# these 2 operators match a values at chebpts1 to the integral function at chebpts1
Il_op = TK1 @ Il @ T_inv
Ir_op = TK1 @ Ir @ T_inv

def to_cheb(phi, l, r):
    phieval = np.array([phi(
        ((r-l)*taui + (r+l))/ 2 
    ) for taui in tau])
    return phieval

def integrate(func, l, r):
    func_v = to_cheb(func, l, r)
    return float(cheb_quad @ func_v) * (r-l) / 2

def integrate_v(vec, l, r):
    return float(cheb_quad @ vec) * (r-l) / 2



def LinearBVPSolver(p, q, f, a, c, zetal0, zetal1, gammal, zetar0, zetar1, gammar,
    C=4.0, TOL=1e-8, iters=40, eval_points=32, force_double=False):

    two_to_C = 2**C

    # solve ui
    Aui = np.array([
        [zetal0, zetal0*a+zetal1, zetal0*a*a + zetal1*2*a],
        [zetar0, zetar0*c+zetar1, zetar0*c*c + zetar1*2*c],
    ])
    Bui = np.array([gammal, gammar])
    uicoeff = solve(Aui, Bui)

    ui  = lambda x: uicoeff[0] + uicoeff[1]*x + uicoeff[2]*x*x
    dui = lambda x: uicoeff[1] + uicoeff[2]*2*x
    ddui= lambda x: uicoeff[2]*2

    ftilde = lambda x: f(x) - (ddui(x) + p(x)*dui(x) + q(x)*ui(x))

    if abs(zetal0) >= abs(zetal1) or abs(zetar0) >= abs(zetar1):
        gl  = lambda x: zetal0 * (x-a) - zetal1
        dgl = lambda x: zetal0
        gr  = lambda x: zetar0 * (x-c) - zetar1
        dgr = lambda x: zetar0
        q0  = 0
    else:
        gl  = lambda x: zetal1 * cosh(x-a) - zetal0 * sinh(x-a)
        dgl = lambda x: zetal1 * sinh(x-a) - zetal0 * cosh(x-a)
        gr  = lambda x: zetar1 * cosh(x-c) - zetar0 * sinh(x-c)
        dgr = lambda x: zetar1 * sinh(x-c) - zetar0 * cosh(x-c)
        q0  = -1

    s = gl(a) * dgr(a) - dgl(a) * gr(a)

    G0 = lambda x, t: ( gl(x)*gr(t)/s if x<=t else gl(t)* gr(x)/s)
    G1 = lambda x, t: (dgl(x)*gr(t)/s if x<=t else gl(t)*dgr(x)/s)
    
    qtilde = lambda x: q(x) - q0

    psil = lambda x: (p(x) * dgr(x) + qtilde(x) * gr(x)) / s
    psir = lambda x: (p(x) * dgl(x) + qtilde(x) * gl(x)) / s

    def PiMatrix(bl, br):
        gl_v = to_cheb(gl, bl, br)
        gr_v = to_cheb(gr, bl, br)
        psil_v = to_cheb(psil, bl, br)
        psir_v = to_cheb(psir, bl, br)
        A = np.zeros((K,K))

        for x in range(K):
            for j in range(K):
                A[x][j] = (
                    psil_v[x] * Il_op[x][j] * gl_v[j] +
                    psir_v[x] * Ir_op[x][j] * gr_v[j] 
                    # see (55) (56)
                )
        A *= (br-bl)/2
        A += np.eye(K)
        return A
    

    class Node:
        def __init__(self, bl, br):
            self.bl = bl
            self.br = br
            self.D = None
            self.E = None
            self.renew()

        def renew(self):
            self.alphal = 0
            self.betal  = 0
            self.deltal = 0
            self.alphar = 0
            self.betar  = 0
            self.deltar = 0
            self.mu     = 1
            self.mul    = 0
            self.mur    = 0

            if not self.leaf():
                self.D.renew()
                self.E.renew()
            
        def leaf(self):
            return self.D == None

        def __len__(self):
            if self.leaf():
                return 1
            else:
                return len(self.D) + len(self.E)

        def split(self):
            if self.leaf():
                midpoint = (self.bl + self.br)/2
                self.D = Node(self.bl, midpoint)
                self.E = Node(midpoint, self.br)
        
        def compute_abd(self):
            nonlocal root
            if self.leaf():
                A = PiMatrix(self.bl, self.br)
                # A_inv = inv(A)
                self.pibar_inv_psil_v   = solve(A, to_cheb(psil, self.bl, self.br))
                self.pibar_inv_psir_v   = solve(A, to_cheb(psir, self.bl, self.br))
                self.pibar_inv_ftilde_v = solve(A, to_cheb(ftilde, self.bl, self.br))

                self.gl_v = to_cheb(gl, self.bl, self.br)
                self.gr_v = to_cheb(gr, self.bl, self.br)
                self.dgl_v = to_cheb(dgl, self.bl, self.br)
                self.dgr_v = to_cheb(dgr, self.bl, self.br)

                gl_v = self.gl_v
                gr_v = self.gr_v

                self.alphal = integrate_v((gl_v * self.pibar_inv_psil_v),   self.bl, self.br)
                self.alphar = integrate_v((gr_v * self.pibar_inv_psil_v),   self.bl, self.br)
                self.betal  = integrate_v((gl_v * self.pibar_inv_psir_v),   self.bl, self.br)
                self.betar  = integrate_v((gr_v * self.pibar_inv_psir_v),   self.bl, self.br)
                self.deltal = integrate_v((gl_v * self.pibar_inv_ftilde_v), self.bl, self.br)
                self.deltar = integrate_v((gr_v * self.pibar_inv_ftilde_v), self.bl, self.br)
                
                
            else:
                self.D.compute_abd()
                self.E.compute_abd()

                Delta = 1 - self.E.alphar * self.D.betal
                if Delta == 0: Delta = 1e-20
                self.alphal = (
                    (1 - self.E.alphal) * (self.D.alphal - self.D.betal * self.E.alphar)
                    / Delta + self.E.alphal
                )
                self.alphar = (
                    self.E.alphar * (1 - self.D.betar) * (1 - self.D.alphal)
                    / Delta + self.D.alphar
                )
                self.betal = (
                    self.D.betal * (1 - self.E.betar) * (1 - self.E.alphal)
                    / Delta + self.E.betal
                )
                self.betar = (
                    (1 - self.D.betar) * (self.E.betar - self.D.betal * self.E.alphar)
                    / Delta + self.D.betar
                )
                # (44) corrected
                self.deltal = (
                    (1 - self.E.alphal) / Delta * self.D.deltal 
                    + self.E.deltal
                    + (self.E.alphal - 1) * self.D.betal / Delta * self.E.deltar 
                )
                self.deltar = (
                    (1 - self.D.betar) / Delta * self.E.deltar
                    + self.D.deltar
                    + (self.D.betar - 1) * self.E.alphar / Delta * self.D.deltal
                )
                # (40) to (45)

        def compute_lambda(self):
            if self.leaf():
                return
            self.D.mu = self.mu
            self.E.mu = self.mu
            self.D.mul = self.mul
            self.E.mur = self.mur
            self.D.mur, self.E.mul = inv(np.array([
                [1, self.E.alphar],
                [self.D.betal, 1],
            ])) @ np.array([
                self.mur * (1 - self.E.betar) - self.mu * self.E.deltar,
                self.mul * (1 - self.D.alphal) - self.mu * self.D.deltal,
            ])
            # (34)
            self.D.compute_lambda()
            self.E.compute_lambda()

        def __call__(self, x, deriv=False):
            if self.leaf():
                return cheb.chebval(
                    2 * (x-self.bl) / (self.br-self.bl) - 1,
                    self.du_coeff if deriv else self.u_coeff
                )
            else:
                midpoint = (self.bl + self.br)/2
                return self.D(x,deriv=deriv) if x < midpoint else self.E(x,deriv=deriv)

                
        def condition(self, S0):
            if ((not self.leaf()) and self.D.leaf() and self.E.leaf()
                and (self.D.Si+self.E.Si <= S0)):
                self.Si = self.D.Si + self.E.Si
                self.D = self.E = None

        def merge_all(self, S0):
            if not self.leaf():
                self.D.merge_all(S0)
                self.E.merge_all(S0)
                self.condition(S0)
        
        def double(self):
            if self.leaf(): self.split()
            else:
                self.D.double()
                self.E.double()

        
    def yield_l(node: Node):
        if node.leaf():
            yield node
        else:
            yield from yield_l(node.D)
            yield from yield_l(node.E)

    def yield_r(node: Node):
        if node.leaf():
            yield node
        else:
            yield from yield_r(node.E)
            yield from yield_r(node.D)

    def compute_J(node: Node):
        Jl = 0
        for Z in yield_l(node):
            Z.Jl = Jl
            Jl += Z.deltal + Z.mul * Z.alphal + Z.mur * Z.betal
        Jr = 0
        for Z in yield_r(node):
            Z.Jr = Jr
            Jr += Z.deltar + Z.mul * Z.alphar + Z.mur * Z.betar

    eval_x = np.random.uniform(a,c, size=eval_points)
    eval_f = np.zeros(eval_points)
    eval_f_last = np.zeros(eval_points)

    root = Node(a, c)

    doubled = False
    for iter_cnt in range(iters+3):
        print(f"iter: {iter_cnt}, len: {len(root)}") 
        root.renew()
        root.compute_abd()
        root.compute_lambda()
        compute_J(root)
        
        # get solution
        Sdiv = 0
        for Z in yield_l(root):
            Z.sigma_v = (
                Z.mu * Z.pibar_inv_ftilde_v 
                + Z.mul * Z.pibar_inv_psil_v
                + Z.mur * Z.pibar_inv_psir_v
            )
            # Z.tau_v = ((Z.br-Z.bl)*tau + (Z.br+Z.bl))/2
            ui_v = to_cheb(ui, Z.bl, Z.br)
            dui_v = to_cheb(dui, Z.bl, Z.br)
            # step 3B
            Pv = (Z.Jl + Il_op @ (Z.gl_v * Z.sigma_v) * (Z.br - Z.bl) / 2) /s
            Qv = (Z.Jr + Ir_op @ (Z.gr_v * Z.sigma_v) * (Z.br - Z.bl) / 2) /s

            Z. u_v =  ui_v + Z. gr_v * Pv + Z. gl_v * Qv
            Z.du_v = dui_v + Z.dgr_v * Pv + Z.dgl_v * Qv
            
            Z.sigmai = T_inv @ Z.sigma_v
            Z. u_coeff = T_inv @ Z. u_v
            Z.du_coeff = T_inv @ Z.du_v
            
            si3, si2, si1 = Z.sigmai[-3:]
            Z.Si = abs(si2) + abs(si1-si3)
            Sdiv = max(Sdiv, Z.Si)

        # now TEST
        for i in range(len(eval_x)):
            # evaluate at random points
            eval_f[i] = root(eval_x[i])

        TEST = np.max(np.abs(eval_f - eval_f_last)) / np.max(np.abs(eval_f + eval_f_last + 1e-30))
        eval_f_last = np.array(eval_f)
        
        if iter_cnt >= iters or len(root) > 1024:
            print(f"doesn't converge after {iter_cnt} steps!!")
            return root

        if TEST < TOL: 
            if doubled or not force_double:
                return root
            else:
                root.double()
                doubled = True
        else:
            doubled = False
            root.merge_all(Sdiv / (2**K))
            for Z in yield_l(root):
                if Z.Si > Sdiv / two_to_C:
                    Z.split()


def NewtonNonlinearSolver(
    f, f_2, f_3, a, c, zetal0, zetal1, gammal, zetar0, zetar1, gammar,
    initial=None, C=4.0, TOL=1e-6, iters=10, eval_points=32, force_double=False,
):

    # solve initial value
    if initial == None:
        Aui = np.array([
            [zetal0, zetal0*a+zetal1, zetal0*a*a + zetal1*2*a],
            [zetar0, zetar0*c+zetar1, zetar0*c*c + zetar1*2*c],
        ])
        Bui = np.array([gammal, gammar])
        uicoeff = solve(Aui, Bui)

        ui  = lambda x: uicoeff[0] + uicoeff[1]*x + uicoeff[2]*x*x
        dui = lambda x: uicoeff[1] + uicoeff[2]*2*x
        ddui= lambda x: uicoeff[2]*2

        u = ui
        du = dui
    else: u, du = initial

    eval_x = np.random.uniform(a,c, size=eval_points)
    eval_f = np.zeros(eval_points)
    eval_f_last = np.zeros(eval_points)

    for cnt in range(iters):
        print(f"Newton iter {cnt}")
        p = lambda x: -f_3(x, u(x), du(x))
        q = lambda x: -f_2(x, u(x), du(x))
        u1_root = LinearBVPSolver(p, q, 
            lambda x: f(x, u(x), du(x)) + p(x) * du(x) + q(x) * u(x),
            a, c, zetal0, zetal1, gammal, zetar0, zetar1, gammar,
            C=C, TOL=TOL, iters=iters*4, 
            eval_points=eval_points, force_double=force_double,
        )

        u = lambda x: u1_root(x)
        du = lambda x: u1_root(x, deriv=True)

        # TEST
        for i in range(len(eval_x)):
            # evaluate at random points
            eval_f[i] = u(eval_x[i])

        TEST = np.max(np.abs(eval_f - eval_f_last)) / np.max(np.abs(eval_f + eval_f_last + 1e-30))
        eval_f_last = np.array(eval_f)

        if TEST < TOL: 
            return u1_root

    print(f"doesn't converge after {iters} steps!!")
    return u1_root


if __name__ == "__main__":

    from scipy.special import jv
    root = LinearBVPSolver(
        lambda x: 1/x,
        lambda x: 1-(100/x)**2,
        lambda x: 0,
        0, 600,
        1,0,0,
        1,0,1,
        TOL=5e-10,
        eval_points=100
    )
    xs = np.linspace(0,600,600)
    plt.plot(xs, [(root(xi) - jv(100,xi)/jv(100,600)) for xi in xs])
    plt.title("Bessel 100 on [0,600]")
    plt.savefig("bessel.jpg")
    plt.show()
    

    root = LinearBVPSolver(
        lambda x: 0,
        lambda x: 1-x**2,
        lambda x: 0,
        -6,6,
        1,0,exp(-18),
        1,0,exp(-18),
        TOL=1e-3,
    )
    xs = np.linspace(-6,6,100)
    plt.plot(xs, [(root(xi) - exp(-xi**2/2)) for xi in xs])
    plt.title("u''+(1-x^2)u=0, exp(-x^2/2)")
    plt.savefig("exp-x2.jpg")
    plt.show()
    
    
    root = LinearBVPSolver(
        lambda x: 0,
        lambda x: 1,
        lambda x: 0,
        0, pi/2,
        0,1,1,
        0,1,0,
    )
    xs = np.linspace(0, pi/2, 100)
    plt.plot(xs, [root(xi)-sin(xi) for xi in xs])
    plt.title("u''+u=0,u'(0)=1,u'(pi/2)=0 (Pure Neumann)")
    plt.savefig("sinx.jpg")
    plt.show()
    

    root = NewtonNonlinearSolver(
        lambda x, u, du: 2*u*du, 
        lambda x, u, du: 2*du, 
        lambda x, u, du: 2*u, 
        -pi/4, pi/4,
        1,0,-1,
        1,0,1,
        initial=(lambda x: x**2, lambda x: 2*x),
        TOL=1e-13
    )
    xs = np.linspace(-pi/4, pi/4, 100)
    plt.plot(xs, [root(xi)-tan(xi) for xi in xs])
    plt.title("u''=2uu', u=tanx")
    plt.savefig("tanx.jpg")
    plt.show()
    
    
    root = NewtonNonlinearSolver(
        lambda x, u, du: du**2 / (2*u), 
        lambda x, u, du: - (du / u)**2/2, 
        lambda x, u, du: du / u, 
        -1, 1,
        1,0,1,
        1,0,1,
        initial = (lambda x: x, lambda x: 1),
        TOL=1e-8,
    )
    xs = np.linspace(-1, 1, 100)
    plt.plot(xs, [root(xi)-xi**2 for xi in xs])
    plt.title("u''=u'^2/(2u), u=x^2")
    plt.savefig("x^2.jpg")
    plt.show()

    root = NewtonNonlinearSolver(
        lambda x, u, du: du * (1-du), 
        lambda x, u, du: 0, 
        lambda x, u, du: 1 - 2 * du, 
        -log(2), log(2),
        0, 1, 1/3,
        1, 0, log(3),
        initial = (lambda x: (x+1)/2, lambda x: 1/2),
        TOL=1e-8,
    )
    xs = np.linspace(-log(2), log(2), 100)
    plt.plot(xs, [root(xi)-log(1+exp(xi)) for xi in xs])
    plt.title("u''=u'(1-u'), u=ln(1+exp(x))")
    plt.savefig("logistic.jpg")
    plt.show()
    
