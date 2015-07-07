r"""
Rijndael-GF

An algebraic implementation of the Rijndael-GF extension of the AES cipher, as
described in [DR02]_. The AES cipher itself is defined to operate on a state
in `(\GF{2})^{8 n_t}` where `n_t \in \{16, 20, 24, 28, 32\}`. Rijndael-GF is a
generalization of AES which allows for operations in `(\GF{2^8})^{n_t}`, 
enabling more algebraically sophisticated study of AES and its variants.
This implementation of Rijndael-GF is suitable for learning purposes, for
comparison to other algebraic ciphers, and for studying various techniques 
of algebraic cryptanalysis of AES. This cipher is different from
:mod:`Mini-AES <sage.crypto.block_cipher.miniaes>`, which is a
teaching tool for beginners to understand the basic structure of AES.

An algebraic implementation of Rijndael-GF is achieved by recognizing that
for every round component function operating on state matrices, `\phi`,
every entry of the output matrix `B = \phi(A)` is representable as a
polynomial with variables being the entries of the input state matrix `A`.
In particular, `B_{r,s} = \sum\limits_i \sum\limits_j \, \alpha_{i,j} A_{i,j},
\, \, \alpha_{i,j} \in \GF{2^8}`, for every entry `B_{r,s}` in `B`.
Correspondingly, this implementation of Rijndael-GF provides a method of the 
form ``phi_poly(i, j)`` for each round component function `\phi` which
returns a polynomial of the above form representing `B_{i,j}`
There are additionally various methods provided which allow for easy
polynomial evaluation and for simple creation of polynomials representing more
complex aspects of the cipher. 

This approach to implementing Rijndael-GF bears some similarity to the 
multivariate quadratic (MQ) systems utilized in :mod:`SR <sage.crypto.mq.sr>`,
in that the MQ systems can be seen as a more specific application of the 
concepts demonstrated in this cipher. Despite this similarity, Rijndael-GF
and :mod:`SR <sage.crypto.mq.sr>` are quite different, as this implementation
seeks to provide a fully generalized algebraic representation of components
of the AES cipher while :mod:`SR <sage.crypto.mq.sr>` is a family of 
parameterizable variants of the AES suitable as a framework for comparing
different cryptanalytic techniques that can be brought to bear on the AES.

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

Each of the round component functions (``shift_rows``, ``mix_columns``, and 
``sub_bytes``) has a ``_poly`` method which takes an index ``i,j`` and returns
a polynomial representation of `B_{i,j}` in terms of entries of `A`. ::

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

    sage: matrix(4, rgf.block_length(), rgf.state_PR.gens())
    [a00 a01 a02 a03]
    [a10 a11 a12 a13]
    [a20 a21 a22 a23]
    [a30 a31 a32 a33]

We can also change the name of the variables by specifying it when we create 
our RijndaelGF object. ::

    sage: rgf = RijndaelGF(4, 6, state_chr='myChr')
    sage: matrix(4, rgf.block_length(), rgf.state_PR.gens())
    [myChr00 myChr01 myChr02 myChr03]
    [myChr10 myChr11 myChr12 myChr13]
    [myChr20 myChr21 myChr22 myChr23]
    [myChr30 myChr31 myChr32 myChr33]

We can evaluate each of these polynomials for a particular input state (in
essence, calculate `B_{i,j}`) as such: ::

    sage: rgf = RijndaelGF(4, 6)
    sage: state = rgf.hex_to_GF('fe7b5170fe7c8e93477f7e4bf6b98071')
    sage: poly = rgf.mix_columns_poly(3, 2, algorithm='decrypt')
    sage: poly(state.list())
    x^7 + x^6 + x^5 + x^2 + x

Each of these ``_poly`` functions can be used as input to ``apply_poly``,
which creates a matrix whose `i,j` th entry equals ``_poly(i, j)``. Then 
``apply_poly`` evaluates each polynomial entry by replacing the
variables with elements from an inputted state, returning the matrix `B`. 
This is equivalent to applying the round component function associated with
this ``_poly`` to a state. ::

    sage: state = rgf.hex_to_GF('c4cedcabe694694e4b23bfdd6fb522fa')
    sage: result = rgf.apply_poly(state, rgf.sub_bytes_poly)
    sage: rgf.GF_to_hex(result)
    '1c8b86628e22f92fb32608c1a8d5932d'
    <BLANKLINE>
    sage: result == rgf.sub_bytes(state)
    True

We can build polynomials over multiple round functions by using
the ``compose`` method. The first and second arguments can be methods like
above which take an index ``i,j`` and return a polynomial representing
`B_{i,j}`. If so, the output will be a function which takes an index ``i,j``
and returns a polynomial representing `B_{i,j}`, where `B = g(f(A))`, `f` is
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
make the first argument, `f`, perform its decryption transformation. Setting
this flag does nothing if the second argument is a function, since ``compose``
returns a function which has its own ``algorithm`` flag. ::

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

In addition to building polynomial representations of state matrices, we can
also build polynomial representations of elements of the expanded key. However,
since the key schedule is defined recursively, it is impossible to build 
polynomials for the key schedule in the same manner as we do for normal state
entries. Instead, an entry of the key schedule is defined by its `i,j` 
positioning in the round key and by the round number of that round key. Such an
element is then represented in terms of entries of the original key matrix
provided at the beginning of encryption or decryption. ::

    sage: rgf.expand_key_poly(1, 2, 0)
    k12
    sage: rgf.expand_key_poly(1, 1, 1)
    k15
    sage: rgf.expand_key_poly(1, 2, 1)
    (x^2 + 1)*k25^254 + (x^3 + 1)*k25^253 + (x^7 + x^6 + x^5 + x^4 + x^3 + 1)*k25^251 + (x^5 + x^2 + 1)*k25^247 + (x^7 + x^6 + x^5 + x^4 + x^2)*k25^239 + k25^223 + (x^7 + x^5 + x^4 + x^2 + 1)*k25^191 + (x^7 + x^3 + x^2 + x + 1)*k25^127 + k10 + (x^6 + x^5 + x)
    
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

    # I'm going to need a series of examples at the top demonstrating how
    # the polynomial representation works and how to use it.
    
    # NOTE: I still need to decide how to implement the key addition.
    # For now, I'll omit it from the algebraic representation.
    def __init__(self, Nb, Nk, state_chr='a', key_chr='k'):
        r"""
        An algebraically generalized version of the AES cipher.
        
        INPUT:
        
        - ``Nb`` -- The block length for this particular instantiation.
        
        - ``Nk`` -- The key length for this particular instantion.

        - ``state_chr`` -- The variable name for polynomials representing
          elements from state matrices.

        - ``key_chr`` -- The variable name for polynomials representing 
          elements of the key schedule.

        EXAMPLES::

            
            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(6, 8)
            sage: rgf
            Rijndael-GF block cipher with block length 6, key length 8, and 14 rounds.
            
        By changing ``state_chr`` we can alter the names of variables in 
        polynomials representing elements from state matrices. ::

            sage: rgf = RijndaelGF(4, 6, state_chr='myChr')
            sage: rgf.mix_columns_poly(3, 2)
            (x + 1)*myChr02 + myChr12 + myChr22 + (x)*myChr32
            
        """
        if Nb not in range(4, 9):
            msg = 'Block length Nb must be in the range 4 - 8, not {0}'
            raise ValueError(msg.format(Nb))
        if Nk not in range(4, 9):
            msg = 'Key length Nk must be in the range 4 - 8, not {0}'
            raise ValueError(msg.format(Nk))

        self._Nb = Nb
        self._Nk = Nk
        self._bits_to_word = 8
        
        from sage.rings.polynomial.polynomial_ring import polygen
        from sage.rings.finite_rings.integer_mod_ring import Integers
        pgen = polygen(Integers(2))
        mod = pgen**8 + pgen**4 + pgen**3 + pgen + 1
        self._F = FiniteField(2**self._bits_to_word, 'x', modulus=mod)
        round_num_table = matrix([[10,11,12,13,14], [11,11,12,13,14],
                                  [12,12,12,13,14], [13,13,13,13,14],
                                  [14,14,14,14,14]])
        self._Nr = round_num_table[self._Nb - 4][self._Nk - 4]

        state_names = [state_chr + str(i) + str(j) 
                       for i in range(4) for j in range(self._Nb)]
        self.state_PR = PolynomialRing(self._F, len(state_names), state_names)
        key_names = [key_chr + str(i) + str(j)
                     for i in range(4) for j in range(self._Nk)]
        ks_names = state_names + key_names
        self.key_PR = PolynomialRing(self._F, len(ks_names), ks_names)
        self._state_vrs = matrix(4, self._Nb, self.state_PR.gens())
        key_vrs_gens = self.key_PR.gens()[len(self.state_PR.gens()):]
        self._key_vrs = matrix(4, self._Nk, key_vrs_gens)

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
        mixcols_E_row = [self._F('x'), self._F('x+1'), self._F('1'), 
                         self._F('1')]
        self._mixcols_E = matrix([mixcols_E_row[-i:] + mixcols_E_row[:-i]
                                  for i in range(4)])
        mixcols_D_row = [self._F('x^3 + x^2 + x'), self._F('x^3 + x + 1'),
                         self._F('x^3 + x^2 + 1'), self._F('x^3 + 1')]
        self._mixcols_D = matrix([mixcols_D_row[-i:] + mixcols_D_row[:-i]
                                  for i in range(4)])

    def __call__(self, text, key, algorithm='encrypt', format='hex'):
        r"""
        Returns a plaintext or ciphertext encryption/decryption of ``text``
        with key ``key``.

        INPUT:

        - ``text`` -- A plaintext to encrypt or a ciphertext to decrypt.

        - ``key`` -- The key to encrypt/decrypt ``text`` with.

        - ``algorithm`` -- Whether to encrypt or decrypt ``text``. Flag for
          encryption is "encrypt", flag for decryption is "decrypt".

        - ``format`` -- The format of ``text`` and ``key``, either "hex" or 
          "binary"

        OUTPUT:

        - The encrypted or decrypted message ``text`` with key ``key``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: text = 'ef053f7c8b3d32fd4d2a64ad3c93071a'
            sage: key = '2d7e86a339d9393ee6570a1101904e16'
            sage: rgf(text, key)
            '84e75b142c8fd5a445312c0a9b2d6699'
            sage: rgf(text, key, algorithm='decrypt')
            '9bf83275406304f050c826ca72d035e6'

        We can use binary strings for ``text`` and ``key``. ::

            sage: text = '11011100011010000011101111011011' * 4
            sage: key = '01000000000011000101101011011110' * 4
            sage: rgf(text, key, format='binary')
            '00011000010110010011100100010111010101001000010010100110101010101111001001100000011111011100100011010001010100110011000111110011'
            sage: rgf(text, key, algorithm='decrypt', format='binary')
            '11000110011001001110000101011101001001010101110001110010000111110000010111111101000011010101101011111100100001010010111000011010'
        """
        
        if algorithm == 'encrypt':
            return self.encrypt(text, key, format)
        elif algorithm == 'decrypt':
            return self.decrypt(text, key, format)
        else:
            raise ValueError(("keyword 'algorithm' must be either 'encrypt' "
                             "or 'decrypt'"))

    def __repr__(self):
        r"""
        Return the string representation of ``self``.

        EXAMPLES ::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(5, 8)
            sage: rgf
            Rijndael-GF block cipher with block length 5, key length 8, and 14 rounds.
        """

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
            raise TypeError("keyword 'H' must be a hex string")

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
        if not isinstance(GF, Matrix_dense) and \
           not isinstance(GF, list) and \
           not is_FiniteFieldElement(GF):
            msg = ("keyword 'GF' must be a matrix over {0}, a list of "
                   "elements from {0}, or a single element from {0}")
            raise TypeError(msg.format(self._F))

        if isinstance(GF, Matrix_dense):
            if not GF.base_ring().characteristic() == 2 or \
               not GF.base_ring().order() == 2**8 or \
               not GF.base_ring().is_field():
                msg = "The elements of keyword 'GF' must all be from"
                raise TypeError(msg.format(self._F))
            return ''.join([self.GF_to_hex(el)
                            for col in GF.columns() for el in col])
        elif isinstance(GF, list):
            if not all([g.parent().characteristic() == 2 and 
                        g.parent().order() == 2**8 and g.parent().is_field()
                        for g in GF]):
                msg = "The elements of keyword 'GF' must all be from"
                raise TypeError(msg.format(self._F))
            return ''.join([self.GF_to_hex(el) for el in GF])
        else:
            if not GF.parent().characteristic() == 2 or \
               not GF.parent().order() == 2**8 or not GF.parent().is_field():
                msg = "keyword  'GF' must be in"
                raise TypeError(msg.format(self._F))
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
            raise TypeError("keyword 'B' must be a binary string")
        
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
            msg = ("keyword 'GF' must be a matrix over {0}, a list of "
                   "elements from {0}, or a single element from {0}")
            raise TypeError(msg.format(self))

        if isinstance(GF, Matrix_dense):
            if not GF.base_ring().characteristic() == 2 or \
               not GF.base_ring().order() == 2**8 or \
               not GF.base_ring().is_field():
                msg = "The elements of keyword 'GF' must all be from"
                raise TypeError(msg.format(self._F))
            return ''.join([self.GF_to_bin(el)
                            for col in GF.columns() for el in col])
        elif isinstance(GF, list):
            if not all([g.parent().characteristic() == 2 and
                        g.parent().order() == 2**8 and g.parent().is_field()
                        for g in GF]):
                msg = "The elements of keyword 'GF' must all be from"
                raise TypeError(msg.format(self._F))
            return ''.join([self.GF_to_bin(el) for el in GF])
        else:
            if not GF.parent().characteristic() == 2 or \
               not GF.parent().order() == 2**8 or not GF.parent().is_field():
                msg = "keyword 'GF' must be in"
                raise TypeError(msg.format(self._F))
            return bin(GF.integer_representation())[2:].zfill(8)

    def encrypt(self, plain, key, format='hex'):
        r"""
        Returns the plaintext ``plain`` encrypted with key ``key`` in the 
        format specified by ``format``.

        INPUT:

        - ``plain`` -- The plaintext to be encrypted.

        - ``key`` -- The key to encrypt ``plain`` with.

        - ``format`` -- The string format of ``key`` and ``plain``, either
          "hex" or "binary".

        OUTPUT:

        - A string of the plaintext ``plain`` encrypted with the key ``key``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: key = 'c81677bc9b7ac93b25027992b0261996'
            sage: plain = 'fde3bad205e5d0d73547964ef1fe37f1'
            sage: rgf.encrypt(plain, key)
            'e767290ddfc6414e3c50a444bec081f0'
            
        We can encrypt binary strings as well. ::

            sage: key = '10010111110000011111011011010001' * 4
            sage: plain = '00000000101000000000000001111011' * 4
            sage: rgf.encrypt(plain, key, format='binary')
            '11010111100100001010001011110010111111001100000001111110010001101110010100000000100011100001000100111011011001000111101111110100'
        """
        if format == 'hex':
            if not isinstance(plain, basestring) or \
               any([c not in '0123456789abcdefABCDEF' for c in plain]):
                raise TypeError("'plain' keyword must be a hex string")
            if len(plain) != 8 * self._Nb:
                msg = "'plain' keyword\'s length must be {0}, not{1}"
                raise ValueError(msg.format(8 * self._Nb, len(plain)))
            if not isinstance(key, basestring) or \
               any([c not in '0123456789abcdefABCDEF' for c in key]):
                raise TypeError("'key' keyword must be a hex string")
            if len(key) != 8 * self._Nk:
                msg = "'key' keyword's length must be {0}, not {1}"
                raise ValueError(msg.format(8 * self._Nk, len(key)))

            state = self.hex_to_GF(plain)
            key_state = self.hex_to_GF(key)
            roundKeys = self.expand_key(key_state)
        elif format == 'binary':
            if not isinstance(plain, basestring) or \
               any([c not in '01' for c in plain]):
                raise TypeError("'plain' keyword must be a binary string")
            if len(plain) != 32 * self._Nb:
                msg = "'plain' keyword's length must be {0}, not {1}"
                raise ValueError(msg.format(32 * self._Nb, len(plain)))
            if not isinstance(key, basestring) or \
               any([c not in '01' for c in key]):
                raise TypeError("'key' keyword must be a binary string")
            if len(key) != 32 * self._Nk:
                msg = "'key' keyword's length must be {0}, not {1}"
                raise ValueError(msg.format(32 * self._Nk, len(key)))

            state = self.bin_to_GF(plain)
            key_state = self.bin_to_GF(key)
            roundKeys = self.expand_key(key_state)
        else:
            raise ValueError(("'format' keyword must be either 'hex' or "
                             "'binary'"))

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
        r"""
        Decrypts ``ciphertext`` with the key ``key`` and returns it as a string
        in the format ``format``.

        INPUT:

        - ``ciphertext`` -- The ciphertext to be decrypted.
        
        - ``key`` -- The key to decrypt ``ciphertext`` with.

        - ``format`` -- The string format that both ``ciphertext`` and ``key``
          must be in, either "hex" or "binary".

        OUTPUT:

        - A string in the format ``format`` of ``ciphertext`` decrypted with
          key ``key``.

        EXAMPLES ::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: key = '2dfb02343f6d12dd09337ec75b36e3f0'
            sage: ciphertext = '54d990a16ba09ab596bbf40ea111702f'
            sage: rgf.decrypt(ciphertext, key)
            '1e1d913b7274ad9b5a4ab1a5f9133b93'

        We can also decrypt messages using binary strings. ::

            sage: key = '00011010000011100011000000111101' * 4
            sage: ciphertext = '00110010001110000111110110000001' * 4
            sage: rgf.decrypt(ciphertext, key, format='binary')
            '10111111101001110011110010101010011111110100001011011000011010000000000000000100000001001110110100001111100011010001101101001011'
        """
        if format == 'hex':
            if not isinstance(ciphertext, basestring) or \
               any([c not in '0123456789abcdefABCDEF' for c in ciphertext]):
                raise TypeError("'ciphertext' keyword must be a hex string")
            if len(ciphertext) != 8 * self._Nb:
                msg = "'ciphertext' keyword's length must be {0}, not{1}"
                raise ValueError(msg.format(8 * self._Nb, len(ciphertext)))
            if not isinstance(key, basestring) or \
               any([c not in '0123456789abcdefABCDEF' for c in key]):
                raise TypeError("'key' keyword must be a hex string")
            if len(key) != 8 * self._Nk:
                msg = "'key' keyword's length must be {0}, not {1}"
                raise ValueError(msg.format(8 * self._Nk, len(key)))

            state = self.hex_to_GF(ciphertext)
            key_state = self.hex_to_GF(key)
            roundKeys = self.expand_key(key_state)
        elif format == 'binary':
            if not isinstance(ciphertext, basestring) or \
               any([c not in '01' for c in ciphertext]):
                raise TypeError(("'ciphertext' keyword must be a binary "
                                 "string"))
            if len(ciphertext) != 32 * self._Nb:
                msg = "'ciphertext' keyword's length must be {0}, not {1}"
                raise ValueError(msg.format(32 * self._Nb, len(ciphertext)))
            if not isinstance(key, basestring) or \
               any([c not in '01' for c in key]):
                raise TypeError("'key' keyword must be a binary string")
            if len(key) != 32 * self._Nk:
                msg = "'key' keyword\'s length must be {0}, not {1}"
                raise ValueError(msg.format(32 * self._Nk, len(key)))

            state = self.bin_to_GF(ciphertext)
            key_state = self.bin_to_GF(key)
            roundKeys = self.expand_key(key_state)
        else:
            raise ValueError(("'format' keyword must be either \'hex\' or "
                             "'binary'"))

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
        r"""
        Raises an error if ``GFm`` is not a proper state matrix.

        INPUT:

        - ``GFm`` -- If ``GFm`` is a proper state matrix, this method returns
          nothing. Otherwise, this method raises an error.

        - ``keyword`` -- The name of ``GFm`` from where this method was called.
          For example, if called from ``sub_bytes``, ``keyword`` would be 
          "state".
          
        - ``expected_cols`` -- The number of columns ``GFm`` should have. 
          Almost always equal to `N_k`, with the exception of ``expand_key``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: bad_state = "I am not a proper state"
            sage: rgf._check_valid_GFmatrix(bad_state, keyword='keyword')
            Traceback (most recent call last):
            ...
            TypeError: keyword 'keyword' must be a 4 x 4 matrix over GF(256)
            <BLANKLINE>
            sage: not_right_size = rgf.hex_to_GF('6353e08c0')
            sage: rgf._check_valid_GFmatrix(not_right_size, keyword='kw')
            Traceback (most recent call last):
            ...
            TypeError: keyword 'kw' must be a 4 x 4 matrix over GF(256)
          
        """
        if expected_cols == None:
            expected_cols = self._Nb
        msg = "keyword '{0}' must be a {1} x {2} matrix over GF({3})"
        msg = msg.format(keyword, 4, expected_cols, self._F.order())
        if not isinstance(GFm, Matrix_dense) or \
           not (GFm.base_ring().order() == self._F.order() and \
                GFm.base_ring().is_field()) or\
           not (GFm.ncols() == expected_cols and GFm.nrows() == 4):
            raise TypeError(msg)

    def _test_poly_input_bounds(self, i, j):
        r"""
        Raises an error if ``i`` or ``j`` are not in the correct bounds.

        INPUT: 

        - ``i`` -- The row index.

        - ``j`` -- The column index.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(6, 4)            
            sage: rgf._test_poly_input_bounds(0, 0)
            sage: rgf._test_poly_input_bounds(5, 0)
            Traceback (most recent call last): 
            ...
            ValueError: keyword 'i' must be between 0 and 4
            sage: rgf._test_poly_input_bounds(0, 7)
            Traceback (most recent call last):
            ...
            ValueError: keyword 'j' must be between 0 and 6
        """
        if i not in range(0, 4):
            raise ValueError("keyword 'i' must be between 0 and 4")
        if j not in range(0, self._Nb):
            msg = "keyword 'j' must be between 0 and {0}"
            raise ValueError(msg.format(self._Nb))

    def expand_key(self, key):
        r"""
        Returns the expanded key schedule from ``key``.

        INPUT:

        - ``key`` -- The key to build a key schedule from. Must be a state 
          matrix of dimensions `4 \times N_k`.

        OUTPUT:

        - A length `Nr` array of `4 \times N_b` matrices corresponding to the
          expanded key. The `n` th entry of the array corresponds to the matrix
          used in the ``add_round_key`` step of the `n` th round.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 6)
            sage: key = '331D0084B176C3FB59CAA0EDA271B565BB5D9A2D1E4B2892'
            sage: key_state = rgf.hex_to_GF(key)
            sage: key_schedule = rgf.expand_key(key_state)
            sage: rgf.GF_to_hex(key_schedule[0])
            '331d0084b176c3fb59caa0eda271b565'
            sage: rgf.GF_to_hex(key_schedule[6])
            '5c5d51c4121f018d0f4f3e408ae9f78c'
        """
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

    def expand_key_poly(self, i, j, round):
        r"""
        Returns a polynomial representing the `i,j` th entry of the ``round``
        th round key.

        INPUT: 

        - ``i`` -- The row position of the element represented by this 
          polynomial.
          
        - ``j`` -- The column position of the element represented by this
          polynomial.

        OUTPUT:

        - A polynomial representing the `i,j` th entry of the ``round`` th 
          round key in terms of entries of the input key. Note that each round
          key has dimensions `4 \times N_b`.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: rgf.expand_key_poly(1, 2, 0)
            k12
            sage: rgf.expand_key_poly(1, 2, 1)
            (x^2 + 1)*k23^254 + (x^3 + 1)*k23^253 + (x^7 + x^6 + x^5 + x^4 + x^3 + 1)*k23^251 + (x^5 + x^2 + 1)*k23^247 + (x^7 + x^6 + x^5 + x^4 + x^2)*k23^239 + k23^223 + (x^7 + x^5 + x^4 + x^2 + 1)*k23^191 + (x^7 + x^3 + x^2 + x + 1)*k23^127 + k10 + k11 + k12 + (x^6 + x^5 + x)

        It should be noted that ``expand_key_poly`` cannot be used with 
        ``apply_poly`` or ``compose``, due to the state polynomials and the
        expanded key polynomials being represented differently.

            sage: rgf.compose(rgf.sub_bytes_poly, rgf.expand_key_poly)
            Traceback (most recent call last):
            ...
            ValueError: expand_key_poly cannot be used with compose
            <BLANKLINE>
            sage: state = rgf.hex_to_GF('00000000000000000000000000000000')
            sage: rgf.apply_poly(state, rgf.expand_key_poly)
            Traceback (most recent call last):
            ...
            ValueError: expand_key_poly cannot be used with apply_poly
        """
        col = round * self._Nb + j
        if col < self._Nk:
            return self._key_vrs[i, col]
        else:

            if col % self._Nk == 0 or (self._Nk > 6 and col % self._Nk == 4):
                # Apply non-linear transformation to col - 1
                recur_r = int((col- 1)/self._Nb)
                recur_j = (col - 1) - (recur_r * self._Nb)
                non_linear = self.expand_key_poly((i+1) % 4, recur_j, recur_r)
                non_linear = self._srd(non_linear)
                non_linear += self._F.gen() ** (int(col / self._Nk) - 1)
                recur_r = int((col- self._Nk)/self._Nb)
                recur_j = (col - self._Nk) - (recur_r * self._Nb)
                return self.expand_key_poly(i, recur_j, recur_r) + non_linear
            else:
                recur_r = int((col - self._Nk)/self._Nb)
                recur_j = (col - self._Nk) - (recur_r * self._Nb)
                result = self.expand_key_poly(i, recur_j, recur_r)
                recur_r = int((col- 1)/self._Nb)
                recur_j = (col - 1) - (recur_r * self._Nb)
                return result + self.expand_key_poly(i, recur_j, recur_r)
        
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
          encryption or decryption. The encryption flag is "encrypt" and the
          decrypt flag is "decrypt".

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
        self._check_valid_GFmatrix(state, 'state')

        if poly_method == self.expand_key_poly:
            raise ValueError("expand_key_poly cannot be used with apply_poly")
        output = matrix(self._F, 4, self._Nb)
        for i in range(4):
            for j in range(4):
                # this is to combat a major performance issue caused by 
                # subbytes
                # Is there a better alternative?
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
          nothing if ``g`` is a function. The encryption flag is "encrypt"
          and the decryption flag is "decrypt".

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
        if f == self.expand_key_poly or g == self.expand_key_poly:
            raise ValueError("expand_key_poly cannot be used with compose")
        if g in self.state_PR:
            f_vals = [f(i, j, algorithm)
                      for i in range(4) for j in range(self._Nb)]
            return g(f_vals)
        else:
            lm = lambda i, j, alg='encrypt': self.compose(f, g(i, j, alg), alg)
            return lm

    def add_round_key(self, state, round_key):
        r"""
        Returns a state matrix over `\GF{2^8}` representing the round key
        addition between state matrices ``state`` and ``round_key``.

        INPUT:
        
        - ``state`` -- The state matrix to have the round key added to.

        - ``round_key`` -- The round key to add to ``state``.

        OUTPUT:

        - A state matrix which is the round key addition of ``state`` and 
          ``round_key``. This transformation is simply the entry-wise addition
          of these two matrices.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: state = rgf.hex_to_GF('36339d50f9b539269f2c092dc4406d23')
            sage: key = rgf.hex_to_GF('7CC78D0E22754E667E24573F454A6531')
            sage: key_schedule = rgf.expand_key(key)
            sage: result = rgf.add_round_key(state, key_schedule[0])
            sage: rgf.GF_to_hex(result)
            '4af4105edbc07740e1085e12810a0812'
        """
        self._check_valid_GFmatrix(state, 'state')
        self._check_valid_GFmatrix(round_key, 'round_key')
        return state + round_key

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
          or as a decryption. The encrypt flag is "encrypt" and the decrypt
          flag is "decrypt".

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
        self._test_poly_input_bounds(i, j)
        if algorithm == 'encrypt':
            var = self._state_vrs[i,j]
            cs = self._sb_E_coeffs
            if no_inversion:
                return sum([cs[i] * (var**(2**i)) for i in range(8)]) + \
                    self._F("x^6 + x^5 + x + 1")
            else:
                return sum([cs[i] * (var**(255 - 2**i)) for i in range(8)]) + \
                    self._F("x^6 + x^5 + x + 1")
        elif algorithm == 'decrypt':
            var = self._state_vrs[i,j]
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
          or the decryption transformation. The encryption flag is "encrypt"
          and the decryption flag is "decrypt".

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
        r"""
        Returns a state matrix in `\GF{2^8}` which is the result of applying
        the SubBytes round component function to ``state``.

        INPUT:

        - ``state`` -- The state matrix to apply SubBytes to.
        
        - ``algorithm`` -- Whether to apply the encryption step of SubBytes
          or its decryption inverse. The flag for encryption is "encrypt" and
          "decrypt" for decryption.

        OUTPUT:

        - The state matrix where SubBytes has been applied to every entry of
          ``state``. The SubBytes transformation consists of first an inversion
          in `\GF{2^8}`, then applying an affine transformation to the result.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: state = rgf.hex_to_GF('d1c4941f7955f40fb46f6c0ad68730ad')
            sage: result = rgf.sub_bytes(state)
            sage: rgf.GF_to_hex(result)
            '3e1c22c0b6fcbf768da85067f6170495'
            <BLANKLINE>
            sage: decryption = rgf.sub_bytes(result, algorithm='decrypt')
            sage: decryption == state
            True
        """
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
          the decryption transformation. The encryption flag is "encrypt" and
          the decryption flag is "decrypt".

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
        self._test_poly_input_bounds(i, j)
        if algorithm == 'encrypt':
            mCols = self._mixcols_E
        elif algorithm == 'decrypt':
            mCols = self._mixcols_D
        else:
            raise ValueError(("keyword 'algorithm' must be either 'encrypt' "
                             "or 'decrypt'"))
        return sum([mCols[i,k] * self._state_vrs[k,j] for k in range(4)])

    def mix_columns(self, state, algorithm='encrypt'):
        r"""
        Returns a state matrix over `\GF{2^8}` to which the round component
        function MixColumns has been applied to.

        INPUT:

        - ``state`` -- The state matrix to apply MixColumns to.
        
        - ``algorithm`` -- Whether to perform the encryption version of 
          MixColumns, or its decryption inverse. The flag for encryption
          is "encrypt" and "decrypt" for decryption.

        OUTPUT:
        
        - A state matrix which is the result of applying MixColumns to 
          ``state``. The MixColumns transformation involves representing each 
          column of ``state`` as a degree 3 polynomial `a(x) \in \GF{2^8}[x]`,
          then replacing that column with another column corresponding to the
          coefficients of `a(x) * (03 x^3 + 01 x^2 + 01 x + 02) (mod x^4 + 1)`.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: state = rgf.hex_to_GF('cd54c7283864c0c55d4c727e90c9a465')
            sage: result = rgf.mix_columns(state)
            sage: rgf.GF_to_hex(result)
            '921f748fd96e937d622d7725ba8ba50c'
            <BLANKLINE>
            sage: decryption = rgf.mix_columns(result, algorithm='decrypt')
            sage: decryption == state
            True
        """
        self._check_valid_GFmatrix(state, 'state')
        return self.apply_poly(state, self.mix_columns_poly, algorithm)

    def shift_rows_poly(self, i, j, algorithm='encrypt'):
        r"""
        Returns a polynomial representing the `i,j` th entry of a state matrix
        after ShiftRows has been applied to it.

        INPUT: 

        - ``i`` -- The row number of the entry represented by this method's
          output.

        - ``j`` -- The column number of the entry represented by this method's 
          output.

        - ``algorithm`` -- Whether to perform ShiftRows encryption step or
          its decryption inverse. The flag for encryption is "encrypt" and
          is "decrypt" for decryption.

        OUTPUT:
        
        - A polynomial in terms of entries of the input state matrix which 
          represents the `i,j` th entry of the output matrix after ShiftRows
          has been applied to it.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: rgf.shift_rows_poly(2, 3)
            a21
            
        We can use this to calculate the `i,j` th entry of a state matrix after
        the decryption version of ShiftRows has been applied to it as such: ::

            sage: poly = rgf.shift_rows_poly(2, 3, algorithm='decrypt')
            sage: state = rgf.hex_to_GF('78c4f708318d3cd69655b701bfc093cf')
            sage: result = poly(state.list())
            sage: rgf.GF_to_hex(result)
            '3c'
            <BLANKLINE>
            sage: output = rgf.shift_rows(state, algorithm='decrypt')
            sage: output[2,3] == result
            True        
        """
        self._test_poly_input_bounds(i, j)
        if algorithm == 'encrypt':
            offs = self._shiftrows_offsets_E
        elif algorithm == 'decrypt':
            offs = self._shiftrows_offsets_D
        else:
            raise ValueError(("keyword 'algorithm' must be either 'encrypt' "
                             "or 'decrypt'"))
        return self._state_vrs[i, (j + offs[4 - self._Nb][i]) % 4]
    
    def shift_rows(self, state, algorithm='encrypt'):
        r"""
        Returns the application of the round component function ShiftRows to
        the state matrix ``state``.

        INPUT:
        
        - ``state`` -- A state matrix over `\GF{2^8}` to which ShiftRows is 
          applied to.

        - ``algorithm`` -- Whether to perform the encryption version of 
          ShiftRows or its decryption inverse. The flag for encryption is 
          "encrypt" and is "decrypt" for decryption.

        OUTPUT:
        
        - A state matrix over `\GF{2^8}` which is the application of ShiftRows
          to ``state``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: state = rgf.hex_to_GF('adcb0f257e9c63e0bc557e951c15ef01')
            sage: result = rgf.shift_rows(state)
            sage: rgf.GF_to_hex(result)
            'ad9c7e017e55ef25bc150fe01ccb6395'
            <BLANKLINE>
            sage: decryption = rgf.shift_rows(result, algorithm='decrypt')
            sage: decryption == state
            True
        """
        self._check_valid_GFmatrix(state, 'state')
        return self.apply_poly(state, self.shift_rows_poly, algorithm)
