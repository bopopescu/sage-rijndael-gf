##########################################################################
#
#       Copyright (C) 2008 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#
##########################################################################


import sage.modular.hecke.hecke_operator
from sage.rings.arith import is_prime
import heilbronn

class HeckeOperator(sage.modular.hecke.hecke_operator.HeckeOperator):
    def apply_sparse(self, x):
        """
        Return the image of x under self. If x is not in self.domain(),
        raise a TypeError.

        EXAMPLES:
            sage: M = ModularSymbols(17,4,-1)
            sage: T = M.hecke_operator(4)
            sage: T.apply_sparse(M.0)
            64*[X^2,(1,8)] + 24*[X^2,(1,10)] - 9*[X^2,(1,13)] + 37*[X^2,(1,16)]
            sage: [ T.apply_sparse(x) == T.hecke_module_morphism()(x) for x in M.basis() ]
            [True, True, True, True]
            sage: N = ModularSymbols(17,4,1)
            sage: T.apply_sparse(N.0)
            Traceback (most recent call last):
            ...
            TypeError: x (=[X^2,(0,1)]) must be in Modular Symbols space of dimension 4 for Gamma_0(17) of weight 4 with sign -1 over Rational Field
        """
        if x not in self.domain():
            raise TypeError("x (=%s) must be in %s"%(x, self.domain()))

        # old version just to check for correctness
        #return self.hecke_module_morphism()(x)

        p = self.index()
        if is_prime(p):
            H = heilbronn.HeilbronnCremona(p)
        else:
            H = heilbronn.HeilbronnMerel(p)

        M = self.parent().module()
        mod2term = M._mod2term
        syms = M.manin_symbols()
        K = M.base_ring()
        R = M.manin_gens_to_basis()

        W = R.new_matrix(nrows=1, ncols = R.nrows())

        B = M.manin_basis()

        #from sage.all import cputime
        #t = cputime()
        v = x.element()
        for i in v.nonzero_positions():
            for h in H:
                entries = syms.apply(B[i], h)
                for k, w in entries:
                    f, s = mod2term[k]
                    if s:
                        W[0,f] += s*K(w)*v[i]

        #print 'sym', cputime(t)
        #t = cputime()
        ans = M( v.parent()((W * R).row(0)) )
        #print 'mul', cputime(t)
        #print 'density: ', len(W.nonzero_positions())/(W.nrows()*float(W.ncols()))

        return ans

