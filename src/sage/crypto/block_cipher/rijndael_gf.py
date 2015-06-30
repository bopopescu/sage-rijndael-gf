r"""
Rijndael-GF

An algebraic implementation of the Rijndael-GF cipher as described in [DR02]_.
AES on its own operates on a state in `(\GF{2})^{8 n_t}` where 
`n_t \in \{16, 20, 24, 28, 32\}`. Rijndael-GF is a generalization of this
which allows for operations in `(\GF{2^8})^{n_t}`, allowing for a more
algebraic treatment of the cipher. The nature of Rijndael-GF 
allows for more algebraically sophisticated study of AES and its variants,
and is suitable for learning purposes, comparison to other algebraic
ciphers, and algebraic cryptanalysis. This cipher is different from
:mod:`Mini-AES <sage.crypto.block_cipher.miniaes>`, which is a teaching tool
for beginners to understand the basic structure of AES.

An algebraic implementation of Rijndael-GF is achieved by recognizing that
for every round component function `\phi`, every entry of the output matrix
`B = \phi(A)` is representable as a polynomial with variables in the 
entries of `A`. In particular, 
`B_{r,s} = \sum\limits_i \sum\limits_j \, \alpha_{i,j} A_{i,j}, \, \,
\alpha_{i,j} \in \GF{2^8}`, for every entry `B_{r,s}` in `B`.
Correspondingly, every round component function in Rijndael-GF has a
method of the form ``method_poly(i, j)`` which returns a polynomial 
representing `B_{i,j}`. There are methods which allow for easy polynomial
evaluation and simple creation of polynomials representing more complex aspects
of the cipher. These concepts bear some similarity to the multivariate 
quadratic (MQ) systems utilized in :mod:`SR <sage.crypto.mq.sr>`, a family of
parameterizable variants of the AES suitable as a framework for comparing
different cryptanalytic techniques that can be brought to bear on the AES. The
concepts demonstrated in this implementation are different from those in
:mod:`SR <sage.crypto.mq.sr>`, however, in that this implementation provides 
a more generalized algebraic representation of the AES, while the MQ system
is a more specific application of these concepts and is primarily intended 
for use as a cryptanalytical tool.

EXAMPLES

We build Rijndael-GF with a block length of 4 and a key length of 6: ::

    sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
    sage: rgf = RijndaelGF(4, 6)

We can perform encryption and decryption by calling the ``encrypt`` and 
``decrypt`` methods, or by calling the Rijndael-GF object explicitly. Note
that the default input format is a hex string. ::

    sage: plaintext = '00112233445566778899aabbccddeeff'
    sage: key = '000102030405060708090a0b0c0d0e0f1011121314151617'
    sage: rgf.encrypt(plaintext, key)
    'dda97ca4864cdfe06eaf70a0ec0d7191'
    sage: rgf.decrypt('dda97ca4864cdfe06eaf70a0ec0d7191', key)
    '00112233445566778899aabbccddeeff'
    <BLANKLINE>
    sage: rgf(plaintext, key)
    'dda97ca4864cdfe06eaf70a0ec0d7191'
    sage: rgf('dda97ca4864cdfe06eaf70a0ec0d7191', key, algorithm='decrypt')
    '00112233445566778899aabbccddeeff'

We can also use binary strings as input. ::

    sage: plain = '11101011100111110000000111001100' * 4
    sage: key = '01100010111101101000110010111010' * 6
    sage: cipher = rgf(plain, key, format='binary')
    sage: cipher
    '11010011000010011010110001000011101110110100110100110010011011111100011011100111110011100111010011001110110100011100000011111011'
    sage: rgf(cipher, key, algorithm='decrypt', format='binary') == plain
    True

Each of the round components has a ``_poly`` method which takes an index 
``i,j`` and returns a polynomial representation of `B_{i,j}` in terms of 
entries of `A`. ::

    sage: i, j = 2, 3
    sage: rgf.sub_bytes_poly(i, j)
    (x^2 + 1)*a23^254 + (x^3 + 1)*a23^253 + (x^7 + x^6 + x^5 + x^4 + x^3 + 1)*a23^251 + (x^5 + x^2 + 1)*a23^247 + (x^7 + x^6 + x^5 + x^4 + x^2)*a23^239 + a23^223 + (x^7 + x^5 + x^4 + x^2 + 1)*a23^191 + (x^7 + x^3 + x^2 + x + 1)*a23^127 + (x^6 + x^5 + x + 1)
    <BLANKLINE>
    sage: rgf.shift_rows_poly(i, j)
    a21
    <BLANKLINE>
    sage: rgf.mix_columns_poly(i, j)
    a03 + a13 + (x)*a23 + (x + 1)*a33

We can see how these variables are organized in `A`: ::

    sage: matrix(4, rgf.block_length(), rgf.P.gens())
    [a00 a01 a02 a03]
    [a10 a11 a12 a13]
    [a20 a21 a22 a23]
    [a30 a31 a32 a33]

We can evaluate each of these polynomials for a particular input state as such:
::

    sage: state = rgf.hex_to_GF('fe7b5170fe7c8e93477f7e4bf6b98071')
    sage: poly = rgf.mix_columns_poly(3, 2, algorithm='decrypt')
    sage: poly(state.list())
    x^7 + x^6 + x^5 + x^2 + x

Each of these ``_poly`` functions can be used as input to ``apply_poly``,
which outputs a matrix whose `i,j` th entry equals that ``_poly`` function
called with `i,j``, which is then evaluated by replacing the variables with 
elements from an inputted state. This is equivalent to applying that 
particular function to a state. ::

    sage: state = rgf.hex_to_GF('c4cedcabe694694e4b23bfdd6fb522fa')
    sage: result = rgf.apply_poly(state, rgf.sub_bytes_poly)
    sage: rgf.GF_to_hex(result)
    '1c8b86628e22f92fb32608c1a8d5932d'
    <BLANKLINE>
    sage: result == rgf.sub_bytes(state)
    True

We can build polynomials like above over multiple round functions by using
the ``compose`` method. The first and second arguments can be methods like
above which take an index ``i,j`` and return a polynomial representing
`B_{i,j}`. The output will be a function which takes an index ``i,j`` and
returns a polynomial representing `B_{i,j}`, where `B = g(f(A))`, `f` is
round function corresponding to the first argument, and `g` is the round
function corresponding with the second argument. ::

    sage: fn = rgf.compose(rgf.shift_rows_poly, rgf.mix_columns_poly)
    sage: fn(2, 1)
    a01 + a12 + (x)*a23 + (x + 1)*a30
    <BLANKLINE>
    sage: state = rgf.hex_to_GF('afb73eeb1cd1b85162280f27fb20d585')
    sage: result = rgf.apply_poly(state, fn)
    sage: new_state = rgf.shift_rows(state)
    sage: new_state = rgf.mix_columns(new_state)
    sage: result == new_state
    True
    <BLANKLINE>
    sage: fn = rgf.compose(rgf.mix_columns_poly, rgf.shift_rows_poly)
    sage: result = rgf.apply_poly(state, fn, algorithm='decrypt')
    sage: new_state = rgf.mix_columns(state, algorithm='decrypt')
    sage: new_state = rgf.shift_rows(new_state, algorithm='decrypt')
    sage: new_state == result
    True

Alternatively, we can use ``compose`` to build polynomials across multiple
round functions like above without having to build our own functions. If we
simply make the second argument a polynomial representing the `i,j` th entry 
of the output state `g(A)`, then ``compose`` will output a polynomial 
representing the `i,j` th entry of `g(f(A))`, where `f` is the round function
corresponding to the first argument. ::

    sage: poly = rgf.mix_columns_poly(0, 3)
    sage: poly
    (x)*a03 + (x + 1)*a13 + a23 + a33
    sage: rgf.compose(rgf.sub_bytes_poly, poly)
    (x^3 + x)*a03^254 + (x^3 + x^2 + x + 1)*a13^254 + (x^2 + 1)*a23^254 + (x^2 + 1)*a33^254 + (x^4 + x)*a03^253 + (x^4 + x^3 + x + 1)*a13^253 + (x^3 + 1)*a23^253 + (x^3 + 1)*a33^253 + (x^7 + x^6 + x^5 + x^3 + 1)*a03^251 + (x^4)*a13^251 + (x^7 + x^6 + x^5 + x^4 + x^3 + 1)*a23^251 + (x^7 + x^6 + x^5 + x^4 + x^3 + 1)*a33^251 + (x^6 + x^3 + x)*a03^247 + (x^6 + x^5 + x^3 + x^2 + x + 1)*a13^247 + (x^5 + x^2 + 1)*a23^247 + (x^5 + x^2 + 1)*a33^247 + (x^7 + x^6 + x^5 + x^4 + x + 1)*a03^239 + (x^2 + x + 1)*a13^239 + (x^7 + x^6 + x^5 + x^4 + x^2)*a23^239 + (x^7 + x^6 + x^5 + x^4 + x^2)*a33^239 + (x)*a03^223 + (x + 1)*a13^223 + a23^223 + a33^223 + (x^6 + x^5 + x^4 + 1)*a03^191 + (x^7 + x^6 + x^2)*a13^191 + (x^7 + x^5 + x^4 + x^2 + 1)*a23^191 + (x^7 + x^5 + x^4 + x^2 + 1)*a33^191 + (x^2 + 1)*a03^127 + (x^7 + x^3 + x)*a13^127 + (x^7 + x^3 + x^2 + x + 1)*a23^127 + (x^7 + x^3 + x^2 + x + 1)*a33^127 + (x^6 + x^5 + x + 1)

We can use ``algorithm='decrypt'`` as an argument to ``compose`` in order to
make the first argument perform its decryption transformation. Setting this
flag does nothing if the second argument is a function. ::

    sage: poly = rgf.shift_rows_poly(2, 1)
    sage: rgf.compose(rgf.mix_columns_poly, poly, algorithm='decrypt')
    (x^3 + x^2 + 1)*a03 + (x^3 + 1)*a13 + (x^3 + x^2 + x)*a23 + (x^3 + x + 1)*a33
    <BLANKLINE>
    sage: state = rgf.hex_to_GF('80121e0776fd1d8a8d8c31bc965d1fee')
    sage: with_decrypt = rgf.compose(rgf.sub_bytes_poly, rgf.shift_rows_poly,
    ....: algorithm='decrypt')
    sage: result_wd = rgf.apply_poly(state, with_decrypt)
    sage: no_decrypt = rgf.compose(rgf.sub_bytes_poly, rgf.shift_rows_poly)
    sage: result_nd = rgf.apply_poly(state, no_decrypt)
    sage: result_wd == result_nd
    True
    
AUTHORS:

- Thomas Gagne (2015-06): initial version

REFERENCES:

.. [DR02] Joan Daemen, Vincent Rijmen. The Design of Rijndael. 
  Springer-Verlag Berlin Heidelberg, 2002.

"""

###########################################################################   
# Copyright (c) 2015 Thomas Gagne <thomasgagne100@gmail.com>        
#                                                                              
# This program is free software; you can redistribute it and/or modify         
# it under the terms of the GNU General Public License as published by         
# the Free Software Foundation; either version 2 of the License, or            
# (at your option) any later version.                                          
#                                                                              
# This program is distributed in the hope that it will be useful,              
# but WITHOUT ANY WARRANTY; without even the implied warranty of               
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                
# GNU General Public License for more details.                                 
#                                                                              
# http://www.gnu.org/licenses/                                                 
###########################################################################

from sage.matrix.constructor import matrix
from sage.matrix.constructor import column_matrix
from sage.matrix.matrix_dense import Matrix_dense
from sage.rings.finite_rings.constructor import FiniteField
from sage.rings.integer import Integer
from sage.structure.sage_object import SageObject
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

class RijndaelGF(SageObject):
    r"""



    """

    # I'm going to need a series of examples at the top demonstrating how
    # the polynomial representation works and how to use it.
    
    # NOTE: I still need to decide how to implement the key addition.
    # For now, I'll omit it from the algebraic representation.
    def __init__(self, Nb, Nk, var_name='a'):
        self._Nb = Nb
        self._Nk = Nk
        self._bits_to_word = 8
        
        from sage.rings.polynomial.polynomial_ring import polygen
        from sage.rings.finite_rings.integer_mod_ring import Integers
        pgen = polygen(Integers(2))
        mod = pgen**8 + pgen**4 + pgen**3 + pgen + 1
        self._F = FiniteField(2**self._bits_to_word, 'x', modulus=mod)
        self._round_num_table = matrix([[10,11,12,13,14], [11,11,12,13,14],
                                       [12,12,12,13,14], [13,13,13,13,14],
                                       [14,14,14,14,14]])
        self._Nr = self._round_num_table[self._Nb - 4][self._Nk - 4]

        names = [var_name + str(i) + str(j) 
                 for i in range(4) for j in range(Nb)]
        self.P = PolynomialRing(self._F, len(names), names)
        self._vrs = matrix(4, self._Nb, self.P.gens())

        self._shiftrows_offsets_E = matrix([[0,1,2,3], [0,1,2,3], [0,1,2,3],
                                   [0,1,2,4], [0,1,3,4]])
        self._shiftrows_offsets_D = matrix([[0,-1,-2,-3], [0,-1,-2,-3],
                                           [0,-1,-2,-3], [0,-1,-2,-4],
                                           [0,-1,-3,-4]])

        self._sb_E_coeffs = [self._F("x^2 + 1"), 
                             self._F("x^3 + 1"), 
                             self._F("x^7 + x^6 + x^5 + x^4 + x^3 + 1"), 
                             self._F("x^5 + x^2 + 1"), 
                             self._F("x^7 + x^6 + x^5 + x^4 + x^2"), 
                             self._F("1"),
                             self._F("x^7 + x^5 + x^4 + x^2 + 1"),
                             self._F("x^7 + x^3 + x^2 + x + 1")]
        self._sb_D_coeffs = [self._F("x^2 + 1"),
                             self._F("x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x"),
                             self._F("x^6 + x^5 + x^4 + x^3 + x^2 + x + 1"), 
                             self._F("x^6 + x^4 + x^3 + x"), 
                             self._F("x^6 + x^5 + x^4 + x^3"),
                             self._F("x^6 + x^4 + x^3 + 1"), 
                             self._F("x^7 + x^6 + x^4 + x^3 + x + 1"),
                             self._F("x^6 + x^5 + x^3 + x^2 + x")]
        self._mixcols_E = matrix(4, 4, [self._F('x'), self._F('x+1'),
                                        self._F('1'), self._F('1'),
                                        self._F('1'), self._F('x'),
                                        self._F('x+1'), self._F('1'),
                                        self._F('1'), self._F('1'), 
                                        self._F('x'), self._F('x+1'),
                                        self._F('x+1'), self._F('1'),
                                        self._F('1'), self._F('x')])
        self._mixcols_D = matrix(4, 4, [self._F('x^3 + x^2 + x'), 
                                        self._F('x^3 + x + 1'),
                                        self._F('x^3 + x^2 + 1'),
                                        self._F('x^3 + 1'),
                                        self._F('x^3 + 1'),
                                        self._F('x^3 + x^2 + x'),
                                        self._F('x^3 + x + 1'),
                                        self._F('x^3 + x^2 + 1'),
                                        self._F('x^3 + x^2 + 1'),
                                        self._F('x^3 + 1'),
                                        self._F('x^3 + x^2 + x'),
                                        self._F('x^3 + x + 1'),
                                        self._F('x^3 + x + 1'),
                                        self._F('x^3 + x^2 + 1'),
                                        self._F('x^3 + 1'),
                                        self._F('x^3 + x^2 + x')])

    def __call__(self, text, key, algorithm='encrypt', format='hex'):
        if algorithm == 'encrypt':
            return self.encrypt(text, key, format)
        elif algorithm == 'decrypt':
            return self.decrypt(text, key, format)
        else:
            raise ValueError(("keyword 'algorithm' must be either 'encrypt' "
                             "or 'decrypt'"))

    def __repr__(self):
        msg = ("Rijndael-GF block cipher with block length {0}, key length "
               "{1}, and {2} rounds.")
        return msg.format(self._Nb, self._Nk, self._Nr)

    def block_length(self):
        r"""
        Returns the block length of this instantiation of Rijndael-GF.

        EXAMPLES::
        
            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 6)
            sage: rgf.block_length()
            4
        """
        return self._Nb

    def key_length(self):
        r"""
        Returns the key length of this instantiation of Rijndael-GF.

        EXAMPLES::
        
            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF 
            sage: rgf = RijndaelGF(4, 8)
            sage: rgf.key_length()
            8
        """
        return self._Nk

    def number_rounds(self):
        r"""
        Returns the number of rounds used in this instantiation of Rijndael-GF.
        Fully determined by the block length and key length of the cipher.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF 
            sage: rgf = RijndaelGF(5, 4)
            sage: rgf.number_rounds()
            11
        """
        return self._Nr
        
    def hex_to_GF(self, H, matrix=True):
        r"""
        Returns the list of or the state matrix of elements of `\GF{2^8}` 
        corresponding to the hex string ``H``. Each element in `\GF{2^8}`
        corresponds to a unique 8-bit binary string, which is represented as
        2 hex characters here. In particular, the element 
        `a_7x^7 + a_6x^6 + a_5x^5 + a_4x^4 + a_3x^3 + a_2x^2 + a_1x^1 + a_0` 
        corresponds to the binary string '`a_7a_6a_5a_4a_3a_2a_1a_0`'.

        INPUT:
        
        - ``H`` -- A hex string where every two hex characters correspond to a
          single value in `\GF{2^8}`
          
        - ``matrix`` -- (default: ``True``) Returns a list if ``False``. 
          Returns a state matrix over `\GF{2^8}` if ``True``.

        OUTPUT:

        - A list of or a state matrix of elements of `\GF{2^8}` where each
          element corresponds to the appropriate hex value in ``H``. 

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: state = rgf.hex_to_GF('1147659047cf663b9b0ece8dfc0bf1f0')
            sage: output = rgf.shift_rows(state)
            sage: rgf.GF_to_hex(output)
            '11cfcef0470ef1909b0b653bfc47668d'

        We can make this method output a list by setting ``matrix`` to 
        ``False``. ::

            sage: rgf.hex_to_GF('2f', matrix=False)
            [x^5 + x^3 + x^2 + x + 1]
            sage: rgf.hex_to_GF('1a2b0f', matrix=False)
            [x^4 + x^3 + x, x^5 + x^3 + x + 1, x^3 + x^2 + x + 1]
        """
        if not isinstance(H, basestring) or \
           any([c not in '0123456789abcdefABCDEF' for c in H]):
            raise TypeError("keyword \'H\' must be a hex string")

        def h_to_gf(h):
            return self._F(map(int, bin(int(h, 16))[2:].zfill(8))[::-1])
        hexes = [H[2*i] + H[2*i+1] for i in range(len(H)/2)]
        result = [h_to_gf(h) for h in hexes]
        if matrix:
            return column_matrix(len(result)/4, 4, result)
        else:
            return result

    def GF_to_hex(self, GF):
        r"""
        Returns the hex string representation of ``GF``.  Each element in 
        `\GF{2^8}` corresponds to a unique 8-bit binary string, which is
        represented as 2 hex characters here. In particular, the element
        `a_7x^7 + a_6x^6 + a_5x^5 + a_4x^4 + a_3x^3 + a_2x^2 + a_1x^1 + a_0`
        corresponds to the binary string '`a_7a_6a_5a_4a_3a_2a_1a_0`'.

        INPUT:
        
        - ``GF`` -- Either a state matrix over `\GF{2^8}`, a list of elements
          from `\GF{2^8}`, or a single element from `\GF{2^8}`

        OUTPUT:

        - A hex string representation of ``GF``, where every two characters in
          the string correspond to a single element in `\GF{2^8}`.

        EXAMPLES::
        
            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4,4)
            sage: F.<a> = GF(2^8)
            sage: els = [a^7 + a^5 + a + 1, a^7 + a^6 + a^4 + a^2]
            sage: rgf.GF_to_hex(els)
            'a3d4'
            <BLANKLINE>
            sage: h = '8fb999c973b26839c7f9d89d85c68c72'
            sage: h == rgf.GF_to_hex(rgf.hex_to_GF(h))
            True

        We can use this to get concise output from round functions. ::
        
            sage: plain = rgf.hex_to_GF('72b86c7c0f0d52d3e0d0da104055036b')
            sage: key = rgf.hex_to_GF('93faa123c2903f4743e4dd83431692de')
            sage: output = rgf.add_round_key(plain, key)
            sage: rgf.GF_to_hex(output)
            'e142cd5fcd9d6d94a3340793034391b5'
        """
        from sage.rings.finite_rings.element_base import is_FiniteFieldElement
        # Make sure to test that the base field of GF is isomorphic to rgf._F
        if not isinstance(GF, Matrix_dense) and \
           not isinstance(GF, list) and \
           not is_FiniteFieldElement(GF):
            msg = ("keyword \'GF\' must be a matrix over {0}, a list of "
                   "elements from {0}, or a single element from {0}")
            raise TypeError(msg.format(self))

        if isinstance(GF, Matrix_dense):
            return ''.join([self.GF_to_hex(el)
                            for col in GF.columns() for el in col])
        elif isinstance(GF, list):
            return ''.join([self.GF_to_hex(el) for el in GF])
        else:
            return hex(GF.integer_representation())[2:].zfill(2)

    def bin_to_GF(self, B, matrix=True):
        r"""
        Returns the list of or the state matrix of elements of `\GF{2^8}`
        corresponding to the binary string ``B``. Each element in `\GF{2^8}`
        corresponds to a unique 8-bit binary string. In particular, 
        the element 
        `a_7x^7 + a_6x^6 + a_5x^5 + a_4x^4 + a_3x^3 + a_2x^2 + a_1x^1 + a_0` 
        corresponds to the binary string '`a_7a_6a_5a_4a_3a_2a_1a_0`'.

        INPUT:

        - ``B`` -- A binary string where every 8 bits correspond to a single
          value in `\GF{2^8}`
          
        - ``matrix`` -- (default: ``True``) Returns a list if ``False``. 
          Returns a state matrix over `\GF{2^8}` if ``True``.

        OUTPUT:

        - A list of or a state matrix of elements of `\GF{2^8}` where each
          element corresponds to the appropriate 8-bit binary string in ``B``. 

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)           
            sage: bs = '11101011100111110000000111001100' * 4
            sage: len(bs)
            128
            sage: state = rgf.bin_to_GF(bs)
            sage: output = rgf.sub_bytes(state)
            sage: rgf.GF_to_bin(output)
            '11101001110110110111110001001011111010011101101101111100010010111110100111011011011111000100101111101001110110110111110001001011'
        
        We can make this method output a list by setting ``matrix`` to
        ``False``. ::        

            sage: bs = '01010011'
            sage: rgf.bin_to_GF(bs, matrix=False)
            [x^6 + x^4 + x + 1]
            sage: bs = '000100000111001000110101110001101101011100110101'
            sage: rgf.bin_to_GF(bs, matrix=False)
            [x^4,
             x^6 + x^5 + x^4 + x,
             x^5 + x^4 + x^2 + 1,
             x^7 + x^6 + x^2 + x,
             x^7 + x^6 + x^4 + x^2 + x + 1,
             x^5 + x^4 + x^2 + 1]
        """
        if not isinstance(B, basestring) or \
           any([c not in '01' for c in B]):
            raise TypeError("keyword \'B\' must be a binary string")
        
        def b_to_gf(b):
            return self._F(map(int, b)[::-1])
        bins = [B[8*i : 8*(i+1)] for i in range(len(B) / 8)]
        result = [b_to_gf(b) for b in bins]
        if matrix:
            return column_matrix(len(result)/4, 4, result)
        else:
            return result

    def GF_to_bin(self, GF):
        r"""
        Returns the binary string representation of ``GF``.  Each element in 
        `\GF{2^8}` corresponds to a unique 8-bit binary string. 
        In particular, the element
        `a_7x^7 + a_6x^6 + a_5x^5 + a_4x^4 + a_3x^3 + a_2x^2 + a_1x^1 + a_0`
        corresponds to the binary string '`a_7a_6a_5a_4a_3a_2a_1a_0`'.

        INPUT:
        
        - ``GF`` -- Either a state matrix over `\GF{2^8}`, a list of elements
          from `\GF{2^8}`, or a single element from `\GF{2^8}`

        OUTPUT:

        - A binary string representation of ``GF``, where every eight 
          characters in the string correspond to a single element in
          `\GF{2^8}`.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: F.<a> = GF(2^8)                                             
            sage: els = [a^7 + a^5 + a + 1, a^7 + a^6 + a^4 + a^2]
            sage: rgf.GF_to_bin(els) 
            '1010001111010100'

        We can use this to get clearer output from the round functions. ::

            sage: plain = '11101011100111110000000111001100' * 4
            sage: plain_state = rgf.bin_to_GF(plain)
            sage: key = '00110011100000001111100111010111' * 4
            sage: key_state = rgf.bin_to_GF(key)
            sage: output = rgf.add_round_key(plain_state, key_state)
            sage: rgf.GF_to_bin(output)                                        
            '11011000000111111111100000011011110110000001111111111000000110111101100000011111111110000001101111011000000111111111100000011011'
        """
        from sage.rings.finite_rings.element_base import is_FiniteFieldElement
        if not isinstance(GF, Matrix_dense) and \
           not isinstance(GF, list) and \
           not is_FiniteFieldElement(GF):
            msg = ("keyword \'GF\' must be a matrix over {0}, a list of "
                   "elements from {0}, or a single element from {0}")
            raise TypeError(msg.format(self))

        if isinstance(GF, Matrix_dense):
            return ''.join([self.GF_to_bin(el)
                            for col in GF.columns() for el in col])
        elif isinstance(GF, list):
            return ''.join([self.GF_to_bin(el) for el in GF])
        else:
            return bin(GF.integer_representation())[2:].zfill(8)

    def encrypt(self, plain, key, format='hex'):
        if format == 'hex':
            if not isinstance(plain, basestring) or \
               any([c not in '0123456789abcdefABCDEF' for c in plain]):
                raise TypeError("\'plain\' keyword must be a hex string")
            if len(plain) != 8 * self._Nb:
                msg = '\'plain\' keyword\'s length must be {0}, not{1}'
                raise ValueError(msg.format(8 * self._Nb, len(plain)))
            if not isinstance(key, basestring) or \
               any([c not in '0123456789abcdefABCDEF' for c in key]):
                raise TypeError("\'key\' keyword must be a hex string")
            if len(key) != 8 * self._Nk:
                msg = '\'key\' keyword\'s length must be {0}, not {1}'
                raise ValueError(msg.format(8 * self._Nk, len(key)))

            state = self.hex_to_GF(plain)
            key_state = self.hex_to_GF(key)
            roundKeys = self.expand_key(key_state)
        elif format == 'binary':
            if not isinstance(plain, basestring) or \
               any([c not in '01' for c in plain]):
                raise TypeError("\'plain\' keyword must be a binary string")
            if len(plain) != 32 * self._Nb:
                msg = '\'plain\' keyword\'s length must be {0}, not {1}'
                raise ValueError(msg.format(32 * self._Nb, len(plain)))
            if not isinstance(key, basestring) or \
               any([c not in '01' for c in key]):
                raise TypeError("\'key\' keyword must be a binary string")
            if len(key) != 32 * self._Nk:
                msg = '\'key\' keyword\'s length must be {0}, not {1}'
                raise ValueError(msg.format(32 * self._Nk, len(key)))

            state = self.bin_to_GF(plain)
            key_state = self.bin_to_GF(key)
            roundKeys = self.expand_key(key_state)
        else:
            raise ValueError(("\'format\' keyword must be either \'hex\' or "
                             "\'binary\'"))

        state = self.add_round_key(state, roundKeys[0])
        for r in range(self._Nr-1):
            state = self.sub_bytes(state, algorithm='encrypt')
            state = self.shift_rows(state, algorithm='encrypt')
            state = self.mix_columns(state, algorithm='encrypt')
            state = self.add_round_key(state, roundKeys[r+1])
        state = self.sub_bytes(state, algorithm='encrypt')
        state = self.shift_rows(state, algorithm='encrypt')
        state = self.add_round_key(state, roundKeys[self._Nr])
            
        if format == 'hex':
            return self.GF_to_hex(state)
        else:
            return self.GF_to_bin(state)

    def decrypt(self, ciphertext, key, format='hex'):
        if format == 'hex':
            if not isinstance(ciphertext, basestring) or \
               any([c not in '0123456789abcdefABCDEF' for c in ciphertext]):
                raise TypeError("\'ciphertext\' keyword must be a hex string")
            if len(ciphertext) != 8 * self._Nb:
                msg = '\'ciphertext\' keyword\'s length must be {0}, not{1}'
                raise ValueError(msg.format(8 * self._Nb, len(ciphertext)))
            if not isinstance(key, basestring) or \
               any([c not in '0123456789abcdefABCDEF' for c in key]):
                raise TypeError("\'key\' keyword must be a hex string")
            if len(key) != 8 * self._Nk:
                msg = '\'key\' keyword\'s length must be {0}, not {1}'
                raise ValueError(msg.format(8 * self._Nk, len(key)))

            state = self.hex_to_GF(ciphertext)
            key_state = self.hex_to_GF(key)
            roundKeys = self.expand_key(key_state)
        elif format == 'binary':
            if not isinstance(ciphertext, basestring) or \
               any([c not in '01' for c in ciphertext]):
                raise TypeError(("\'ciphertext\' keyword must be a binary "
                                 "string"))
            if len(ciphertext) != 32 * self._Nb:
                msg = '\'ciphertext\' keyword\'s length must be {0}, not {1}'
                raise ValueError(msg.format(32 * self._Nb, len(ciphertext)))
            if not isinstance(key, basestring) or \
               any([c not in '01' for c in key]):
                raise TypeError("\'key\' keyword must be a binary string")
            if len(key) != 32 * self._Nk:
                msg = '\'key\' keyword\'s length must be {0}, not {1}'
                raise ValueError(msg.format(32 * self._Nk, len(key)))

            state = self.bin_to_GF(ciphertext)
            key_state = self.bin_to_GF(key)
            roundKeys = self.expand_key(key_state)
        else:
            raise ValueError(("\'format\' keyword must be either \'hex\' or "
                             "\'binary\'"))

        state = self.add_round_key(state, roundKeys[self._Nr])
        state = self.shift_rows(state, algorithm='decrypt')
        state = self.sub_bytes(state, algorithm='decrypt')
        for r in range(self._Nr-1):
            state = self.add_round_key(state, roundKeys[self._Nr - r - 1])
            state = self.mix_columns(state, algorithm='decrypt')
            state = self.shift_rows(state, algorithm='decrypt')
            state = self.sub_bytes(state, algorithm='decrypt')
        state = self.add_round_key(state, roundKeys[0])

        if format == 'hex':
            return self.GF_to_hex(state)
        else:
            return self.GF_to_bin(state)
                  
    def _check_valid_GFmatrix(self, GFm, keyword, expected_cols=None):
        if expected_cols == None:
            expected_cols = self._Nb
        msg = "keyword \'{0}\' must be a {1} x {2} matrix over GF({3})"
        msg = msg.format(keyword, 4, expected_cols, self._F.order())
        if not isinstance(GFm, Matrix_dense) or \
           not (GFm.base_ring().order() == self._F.order() and \
                GFm.base_ring().is_field()) or \
           not (GFm.ncols() == expected_cols and GFm.nrows() == 4):
            raise TypeError(msg)

    # I'm going to keep this code here until I'm certain I don't want vectors..
    def _check_valid_GFvector(self, GFv, keyword, expected_length=None):
        if expected_length == None:
            expected_length = self._Nb * 4
        msg = "keyword \'{0}\' must be a length {1} vector over GF({2})"
        msg = msg.format(keyword, expected_length, self._F.order())
        if not all([el.parent().order() == self._F.order() and 
                    el.parent().is_field() for el in GFv]) or \
           not len(GFv) == expected_length:
            raise TypeError(msg)

    # I'm going to need to decide how I want to deal with keys...
    def expand_key(self, key):
        self._check_valid_GFmatrix(key, 'key', self._Nk)
        # Is this bad form?
        def add_cols(col1, col2):
            return map(lambda (x,y): x + y, zip(col1, col2))
        
        key_cols = []
        for i in range(self._Nb * (self._Nr + 1)):
            key_cols.append([])

        for j in range(self._Nk):
            key_cols[j] = list(key.columns()[j])
        for j in range(self._Nk, self._Nb * (self._Nr + 1)):
            if j % self._Nk == 0:
                # Apply non-linear function to k[j - 1]
                add_key = map(self._srd, key_cols[j - 1])
                add_key = add_key[1:] + add_key[:1]
                add_key[0] += self._F.gen() ** (int(j / self._Nk) - 1)
                key_cols[j] = add_cols(key_cols[j - self._Nk], add_key)
            else:
                add_key = key_cols[j - 1]
                if self._Nk > 6 and j % self._Nk == 4:
                    add_key = map(self._srd, add_key)
                key_cols[j] = add_cols(key_cols[j - self._Nk], add_key)

        # Copy the expanded columns into 4xNb blocks, ready for addRoundKey
        round_keys = []
        for r in range(self._Nr + 1):
            rk = column_matrix([key_cols[r*self._Nb + i] 
                                for i in range(self._Nb)])
            round_keys.append(rk)
        return round_keys

    def apply_poly(self, state, poly_method, algorithm='encrypt'):
        r"""
        Returns the state matrix in `\GF{2^8}` where the `i,j` th entry of the
        output equals the polynomial ``poly_method(i, j, algorithm)`` 
        evaluated by setting its variables equal to the corresponding
        entries of ``state``.

        INPUT:
        
        - ``state`` -- A state matrix over `\GF{2^8}` to which ``poly_method``
          is applied to.

        - ``poly_method`` -- The method to apply to the elements of ``state``.
          ``poly_method`` must be callable as ``poly_method(i, j, algorithm)``
          and must return a polynomial.

        - ``algorithm`` -- Passed directly to ``poly_method`` to determine
          encryption or decryption. The encryption flag is 'encrypt' and the
          decrypt flag is 'decrypt'.

        OUTPUT:

        - A state matrix in `\GF{2^8}` where the `i,j` th element of that 
          matrix equals ``poly_method(i, j, algorithm)``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: state = rgf.hex_to_GF('3b59cb73fcd90ee05774222dc067fb68')
            sage: result = rgf.apply_poly(state, rgf.shift_rows_poly)
            sage: rgf.GF_to_hex(result)
            '3bd92268fc74fb735767cbe0c0590e2d'

        Calling ``apply_poly`` with the poly_method of a round component
        (e.g. ``sub_bytes_poly``) is identical to calling that round 
        component method itself. ::

            sage: state = rgf.hex_to_GF('4915598f55e5d7a0daca94fa1f0a63f7')
            sage: ap_result = rgf.apply_poly(state, rgf.sub_bytes_poly)
            sage: direct_result = rgf.sub_bytes(state)
            sage: direct_result == ap_result
            True
        
        """
        # I should include an example here using link()....
        self._check_valid_GFmatrix(state, 'state')

        output = matrix(self._F, 4, self._Nb)
        for i in range(4):
            for j in range(4):
                # this is to combat a major performance issue caused by 
                # subbytes
                # Is there a better alternative?
                # This might also cause issues with key addition....
                if poly_method == self.sub_bytes_poly:
                    if algorithm == 'encrypt':
                        p = poly_method(i, j, algorithm)
                        p = p(state.list())
                    else:
                        p = poly_method(i, j, algorithm, no_inversion=True)
                        p = p(state.list()) ** 254
                else:
                    p = poly_method(i, j, algorithm)
                    p = p(state.list())
                output[i,j] = p
        return output

    def compose(self, f, g, algorithm='encrypt'):
        r"""
        If ``g`` is a method which takes an index ``i,j`` and returns a 
        polynomial representing the `i,j` th entry of the matrix `g(A)` in
        terms of entries of `A`, where `A` is a state matrix, then 
        ``compose(f,g)`` returns a method which takes an index ``i,j`` and 
        returns a polynomial representing the `i,j` th entry of the matrix
        `g(f(A))` in terms of entries of `A`.
        If ``g`` is a polynomial representing the `i,j` th entry of the matrix
        `g(A)` in terms of entries of `A`, then ``compose(f,g)`` returns a 
        polynomial representing the `i,j` th entry of `g(f(A))` in terms of
        entries of `A`.

        INPUT:

        - ``f`` -- A function of the form ``f(i, j, algorithm)`` which returns
          a polynomial reprentation of the `i,j` th entry of `f(A)`.

        - ``g`` -- A function of the form ``g(i, j, algorithm)`` which returns
          a polynomial reprentation of the `i,j` th entry of `g(A)` OR a 
          polynomial representing the `i,j` th entry of `g(A)`.

        - ``algorithm`` -- Whether ``f`` and ``g`` should use their
          encryption transformations or their decryption transformations. Does
          nothing if ``g`` is a function. The encryption flag is 'encrypt'
          and the decryption flag is 'decrypt'.

        OUTPUT:

        - If ``g`` is a method, this returns a new method which takes an
          index ``i,j`` and returns a polynomial representing the `i,j` th 
          entry of `g(f(A))`. If ``g`` is a polynomial, this returns a
          polynomial representing the `i,j` th entry of `g(f(A))`.

        EXAMPLES

        This function allows us to determine the polynomial representations
        of entries across multiple round functions. For example, if we
        wanted a polynomial representing the ``1,3`` entry of a matrix after
        we first apply ShiftRows and then MixColumns to that matrix, we do: ::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: mcp = rgf.mix_columns_poly(1, 3); mcp
            a03 + (x)*a13 + (x + 1)*a23 + a33
            sage: result = rgf.compose(rgf.shift_rows_poly, mcp)
            sage: result
            a03 + (x)*a10 + (x + 1)*a21 + a32

        We can test the correctness of this: ::

            sage: state = rgf.hex_to_GF('fa636a2825b339c940668a3157244d17')
            sage: new_state = rgf.shift_rows(state)
            sage: new_state = rgf.mix_columns(new_state)
            sage: result(state.list()) == new_state[1,3]
            True

        We can also use ``compose`` to make a new function which returns the 
        polynomial representing the `i,j` th entry of a matrix after multiple
        round functions: ::

            sage: fn = rgf.compose(rgf.shift_rows_poly, rgf.mix_columns_poly)
            sage: fn(1, 3)
            a03 + (x)*a10 + (x + 1)*a21 + a32 
            <BLANKLINE>
            sage: fn2 = rgf.compose(rgf.sub_bytes_poly, fn)

        If we use ``compose`` to make a new function, we can use that function
        as input to ``apply_poly``. ::
        
            sage: state = rgf.hex_to_GF('36400926f9336d2d9fb59d23c42c3950')
            sage: result = rgf.apply_poly(state, fn)
            sage: rgf.GF_to_hex(result)
            'f4bcd45432e554d075f1d6c51dd03b3c'
            <BLANKLINE>
            sage: new_state = rgf.shift_rows(state)
            sage: new_state = rgf.mix_columns(new_state)
            sage: result == new_state
            True
        """
        if g in self.P:
            f_vals = [f(i, j, algorithm)
                      for i in range(4) for j in range(self._Nb)]
            return g(f_vals)
        else:
            lm = lambda i, j, alg='encrypt': self.compose(f, g(i, j, alg), alg)
            return lm
                
    def add_round_key(self, state, round_key):
        self._check_valid_GFmatrix(state, 'state')
        self._check_valid_GFmatrix(round_key, 'round_key')
        return state + round_key

    # no_inversion is here because the order of computations during actual
    # encrpytion are *really* slowed down if we try to do inversion before
    # substituting in our arguments. 
    def sub_bytes_poly(self, i, j, algorithm='encrypt', no_inversion=False):
        r"""
        Returns a polynomial representing the `i,j` th entry of a state matrix
        after the SubBytes method has been applied to it. 

        INPUT:

        - ``i`` -- The row number of the entry represented by this method's
          output.

        - ``j`` -- The column number of the entry represented by this 
          method's output.

        - ``algorithm`` -- Whether to return the polynomial as an encryption
          or as a decryption. The encrypt flag is 'encrypt' and the decrypt
          flag is 'decrypt'.

        - ``no_inversion`` -- Don't include the inversion step, only perform 
          the affine transformation.

        OUTPUT:

        - A polynomial representing the `i,j` th entry of a state matrix after 
          the SubBytes method has been applied to it.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: rgf.sub_bytes_poly(2, 3)
            (x^2 + 1)*a23^254 + (x^3 + 1)*a23^253 + (x^7 + x^6 + x^5 + x^4 + x^3 + 1)*a23^251 + (x^5 + x^2 + 1)*a23^247 + (x^7 + x^6 + x^5 + x^4 + x^2)*a23^239 + a23^223 + (x^7 + x^5 + x^4 + x^2 + 1)*a23^191 + (x^7 + x^3 + x^2 + x + 1)*a23^127 + (x^6 + x^5 + x + 1)

        We can use this polynomial to calculate the `i,j` th entry of the 
        output matrix for any given state as such: ::

            sage: state = rgf.hex_to_GF('6385b79ffc538df997be478e7547d691')
            sage: poly = rgf.sub_bytes_poly(2, 3)
            sage: poly(state.list())
            x^7 + x^6 + x^5 + x^4 + x^2 + x

        We can set ``no_inversion`` to ``True`` to get a polynomial
        representation of solely the affine transformation. ::

            sage: rgf.sub_bytes_poly(0, 2, no_inversion=True)
            (x^7 + x^3 + x^2 + x + 1)*a02^128 + (x^7 + x^5 + x^4 + x^2 + 1)*a02^64 + a02^32 + (x^7 + x^6 + x^5 + x^4 + x^2)*a02^16 + (x^5 + x^2 + 1)*a02^8 + (x^7 + x^6 + x^5 + x^4 + x^3 + 1)*a02^4 + (x^3 + 1)*a02^2 + (x^2 + 1)*a02 + (x^6 + x^5 + x + 1)

        When generating a decryption polynomial, determining the inverse
        polynomial of the polynomial representing the affine transformation
        results in a polynomial with thousands of terms which can be very slow
        to calculate. In order to speed up decryption, it is recommended to
        calculate the decryption polynomial with ``no_inversion=True``, 
        evaluate the arguments, then perform the inversion AFTER this result
        has been calculated. ::

            sage: poly = rgf.sub_bytes_poly(0, 0, 
            ....: algorithm='decrypt', no_inversion=True)
            sage: state = rgf.hex_to_GF('b415f8016858552e4bb6124c5f998a4c')
            sage: poly(state.list()) ^ -1
            x^7 + x^6 + x^2 + x

        We can easily apply the SubBytes transformation to a single value
        by placing it into the first entry of an otherwise empty matrix. ::

            sage: state = rgf.hex_to_GF('3f000000000000000000000000000000')
            sage: poly = rgf.sub_bytes_poly(0, 0)
            sage: rgf.GF_to_hex(poly(state.list()))
            '75'
        """
        
        #checks here -- make sure to check that no_inversion is a bool, etc.

        if algorithm == 'encrypt':
            var = self._vrs[i,j]
            cs = self._sb_E_coeffs
            if no_inversion:
                return sum([cs[i] * (var**(2**i)) for i in range(8)]) + \
                    self._F("x^6 + x^5 + x + 1")
            else:
                return sum([cs[i] * (var**(255 - 2**i)) for i in range(8)]) + \
                    self._F("x^6 + x^5 + x + 1")
        elif algorithm == 'decrypt':
            var = self._vrs[i,j]
            cs = self._sb_D_coeffs
            result = (sum([cs[i] * var**(2**i) for i in range(8)]) + \
                        self._F("x^2 + 1"))
            if no_inversion:
                return result
            else:
                return result ** 254
        else:
            raise ValueError(("keyword 'algorithm' must be either 'encrypt' "
                             "or 'decrypt'"))

    def _srd(self, el, algorithm='encrypt'):
        r"""
        Returns the application of the non-linear transformation SubBytes
        (also known as SRD) to ``el``.

        INPUT: 

        - ``el`` -- An element of `\GF{2^8}`.

        - ``algorithm`` -- Whether to perform the encryption transformation
          or the decryption transformation. The encryption flag is 'encrypt'
          and the decryption flag is 'decrypt'.

        OUTPUT:
        
        - The result of the application of the non-linear transformation 
          SubBytes to ``el``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: el = rgf.hex_to_GF('2A', matrix=False)[0]
            sage: rgf._srd(el)
            x^7 + x^6 + x^5 + x^2 + 1
        """

        if algorithm == 'encrypt':
            p = self.sub_bytes_poly(0, 0, algorithm)
            return (p([el] + [self._F.zero()]*((4 * self._Nb)-1)))
        elif algorithm == 'decrypt':
            p = self.sub_bytes_poly(0, 0, algorithm, no_inversion=True)
            return (p([el] + [self._F.zero()]*((4 * self._Nb)-1))) ** 254
        else:
            raise ValueError(("keyword 'algorithm' must be either 'encrypt' "
                             "or 'decrypt'"))

    def sub_bytes(self, state, algorithm='encrypt'):
        self._check_valid_GFmatrix(state, 'state')
        return self.apply_poly(state, self.sub_bytes_poly, algorithm)

    def mix_columns_poly(self, i, j, algorithm='encrypt'):
        r"""
        Returns a polynomial representing the `i,j` th entry of a state matrix
        after the MixColumns method has been applied to it.

        INPUT: 
        
        - ``i`` -- The row number of the entry represented by this method's 
          output.
          
        - ``j`` -- The column number of the entry represented by this method's
          output.

        - ``algorithm`` -- Whether to perform the encryption transformation or
          the decryption transformation. The encryption flag is 'encrypt' and
          the decryption flag is 'decrypt'.

        OUTPUT:
        
        - A polynomial in terms of entries of the input state matrix which 
          represents the `i,j` th entry of the output matrix after MixColumns 
          has been applied to it.

        EXAMPLES::
        
            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: rgf.mix_columns_poly(3, 1)
            (x + 1)*a01 + a11 + a21 + (x)*a31

        We can use this to calculate the `i,j` th entry of a state matrix after
        the decryption version of MixColumns has been applied to it as such: ::

            sage: poly = rgf.mix_columns_poly(2, 2, algorithm='decrypt')
            sage: state = rgf.hex_to_GF('a761ca9b97be8b45d8ad1a611fc97369')
            sage: result = poly(state.list())
            sage: rgf.GF_to_hex(result)
            'b7'
            <BLANKLINE>
            sage: output = rgf.mix_columns(state, algorithm='decrypt')
            sage: output[2,2] == result
            True
        """
        if algorithm == 'encrypt':
            mCols = self._mixcols_E
        elif algorithm == 'decrypt':
            mCols = self._mixcols_D
        else:
            raise ValueError(("keyword 'algorithm' must be either 'encrypt' "
                             "or 'decrypt'"))
        return sum([mCols[i,k] * self._vrs[k,j] for k in range(4)])

    def mix_columns(self, state, algorithm='encrypt'):
        self._check_valid_GFmatrix(state, 'state')
        return self.apply_poly(state, self.mix_columns_poly, algorithm)

    def shift_rows_poly(self, i, j, algorithm='encrypt'):
        if algorithm == 'encrypt':
            offs = self._shiftrows_offsets_E
        elif algorithm == 'decrypt':
            offs = self._shiftrows_offsets_D
        else:
            raise ValueError(("keyword 'algorithm' must be either 'encrypt' "
                             "or 'decrypt'"))
        return self._vrs[i, (j + offs[4 - self._Nb][i]) % 4]
    
    def shift_rows(self, state, algorithm='encrypt'):
        self._check_valid_GFmatrix(state, 'state')
        return self.apply_poly(state, self.shift_rows_poly, algorithm)
