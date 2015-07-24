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
for each round component function `\phi` of AES (SubBytes, ShiftRows, etc.) 
operating on state matrices, every entry of the output matrix `B = \phi(A)` is
representable as a polynomial with variables being the entries of the input
state matrix `A`. Correspondingly, this implementation of Rijndael-GF provides
a method of the form ``phi_poly(row, col)`` for each round component function
`\phi` which returns a polynomial representing `B_{i,j} = \phi(A)_{i,j}` in
terms of the entries of `A`. There are additionally various methods provided 
which allow for easy polynomial evaluation and for simple creation of
polynomials representing more complex aspects of the cipher. 

This approach to implementing Rijndael-GF bears some similarity to the 
multivariate quadratic (MQ) systems utilized in :mod:`SR <sage.crypto.mq.sr>`,
in that the MQ systems can be seen as a more specific application of the 
concepts demonstrated in this cipher. Despite this similarity, Rijndael-GF
and :mod:`SR <sage.crypto.mq.sr>` are quite different, as this implementation
seeks to provide a fully generalized algebraic representation of the 
individual components of as well as the whole AES cipher while 
:mod:`SR <sage.crypto.mq.sr>` is a family of parameterizable variants of the
AES suitable as a framework for comparing different cryptanalytic techniques
that can be brought to bear on the AES.

EXAMPLES

We build Rijndael-GF with a block length of 4 and a key length of 6: ::

    sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
    sage: rgf = RijndaelGF(4, 6)

We can encrypt and decrypt plaintexts and ciphertexts by calling the 
``encrypt`` and ``decrypt`` methods or by calling the Rijndael-GF object
explicitly. Note that the default input format is a hex string. ::

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

We can also use binary strings as input and output. ::

    sage: plain = '11101011100111110000000111001100' * 4
    sage: key = '01100010111101101000110010111010' * 6
    sage: ciphertext = rgf(plain, key, format='binary')
    sage: ciphertext
    '11010011000010011010110001000011101110110100110100110010011011111100011011100111110011100111010011001110110100011100000011111011'
    sage: rgf(ciphertext, key, algorithm='decrypt', format='binary') == plain
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

We can see how these polynomial terms are organized in `A`: ::

    sage: matrix(4, rgf.block_length(), rgf._state_PR.gens())
    [a00 a01 a02 a03]
    [a10 a11 a12 a13]
    [a20 a21 a22 a23]
    [a30 a31 a32 a33]

We can evaluate each of these polynomials for a particular input state (in
essence, calculate `B_{i,j}`) as such: ::

    sage: rgf = RijndaelGF(4, 6)
    sage: state = rgf.hex_to_GF('fe7b5170fe7c8e93477f7e4bf6b98071')
    sage: poly = rgf.mix_columns_poly(3, 2, algorithm='decrypt')
    sage: poly(state.list())
    x^7 + x^6 + x^5 + x^2 + x

We can use the ``apply_poly`` method to build a matrix whose `i,j` th
entry equals the polynomial ``_poly(i, j)`` evaluated by replacing each
variable with the corresponding entry from a particular state matrix ``state``.
Essentially, ``apply_poly`` calculates `B = \phi(A)`, where `A` is specified
as ``state`` and `\phi` is the round component function corresponding to the
``_poly`` method passed to ``apply_poly``. Calling ``apply_poly`` is also
equivalent to applying the round component function associated with ``_poly``
to a state. ::

    sage: state = rgf.hex_to_GF('c4cedcabe694694e4b23bfdd6fb522fa')
    sage: result = rgf.apply_poly(state, rgf.sub_bytes_poly)
    sage: rgf.GF_to_hex(result)
    '1c8b86628e22f92fb32608c1a8d5932d'
    <BLANKLINE>
    sage: result == rgf.sub_bytes(state)
    True

Just like the ``_poly`` methods which correspond to a round component
function, we can build our own ``_poly`` methods which correspond to
the composition of multiple round component functions with the ``compose``
method. In this case, the first and second arguments to ``compose`` are
``_poly`` methods which take an index ``i,j`` and return a polynomial
representing `B_{i,j}`. With these arguments, the output will be a function
which takes an index ``i,j`` and returns a polynomial representing `B_{i,j}`, 
where `B = g(f(A))`, where `f` is the round function corresponding to the first
argument, and where `g` is the round function corresponding with the second
argument. ::

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

Alternatively, we can use ``compose`` to build polynomials representing entries
after multiple round component functions like above without having to 
explicitly build our own functions. If we simply make the second argument a
polynomial representing the `i,j` th entry of the output state `g(A)`, then 
``compose`` will output a polynomial representing the `i,j` th entry of
`g(f(A))`, where `f` is the round function corresponding to the first
argument. ::

    sage: poly = rgf.mix_columns_poly(0, 3)
    sage: poly
    (x)*a03 + (x + 1)*a13 + a23 + a33
    sage: rgf.compose(rgf.sub_bytes_poly, poly)
    (x^3 + x)*a03^254 + (x^3 + x^2 + x + 1)*a13^254 + (x^2 + 1)*a23^254 + (x^2 + 1)*a33^254 + (x^4 + x)*a03^253 + (x^4 + x^3 + x + 1)*a13^253 + (x^3 + 1)*a23^253 + (x^3 + 1)*a33^253 + (x^7 + x^6 + x^5 + x^3 + 1)*a03^251 + (x^4)*a13^251 + (x^7 + x^6 + x^5 + x^4 + x^3 + 1)*a23^251 + (x^7 + x^6 + x^5 + x^4 + x^3 + 1)*a33^251 + (x^6 + x^3 + x)*a03^247 + (x^6 + x^5 + x^3 + x^2 + x + 1)*a13^247 + (x^5 + x^2 + 1)*a23^247 + (x^5 + x^2 + 1)*a33^247 + (x^7 + x^6 + x^5 + x^4 + x + 1)*a03^239 + (x^2 + x + 1)*a13^239 + (x^7 + x^6 + x^5 + x^4 + x^2)*a23^239 + (x^7 + x^6 + x^5 + x^4 + x^2)*a33^239 + (x)*a03^223 + (x + 1)*a13^223 + a23^223 + a33^223 + (x^6 + x^5 + x^4 + 1)*a03^191 + (x^7 + x^6 + x^2)*a13^191 + (x^7 + x^5 + x^4 + x^2 + 1)*a23^191 + (x^7 + x^5 + x^4 + x^2 + 1)*a33^191 + (x^2 + 1)*a03^127 + (x^7 + x^3 + x)*a13^127 + (x^7 + x^3 + x^2 + x + 1)*a23^127 + (x^7 + x^3 + x^2 + x + 1)*a33^127 + (x^6 + x^5 + x + 1)

We can use ``algorithm='decrypt'`` as an argument to ``compose`` in order to
make the first argument, `f`, perform its decryption transformation, provided
the second argument is a polynomial. Setting this flag does nothing if the
second argument is a function, since ``compose`` returns a function which must
have the same arguments as the other ``_poly(row, col, algorithm='encrypt')`` 
methods. ::

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

We can see how key variables are organized in the original key below. Note
that because key variables are subkey entries, if the key length is longer
than the block length we will have subkey entries from multiple rounds in
the original key matrix. ::

    sage: rgf.key_vrs
    [k000 k001 k002 k003 k100 k101]
    [k010 k011 k012 k013 k110 k111]
    [k020 k021 k022 k023 k120 k121]
    [k030 k031 k032 k033 k130 k131]

In addition to building polynomial representations of state matrices, we can
also build polynomial representations of elements of the expanded key with the
``expand_key_poly`` method. However, since the key schedule is defined
recursively, it is impossible to build polynomials for the key schedule in
the same manner as we do for normal state entries. Instead, ``expand_key_poly``
takes an index ``i,j`` and a round number ``round``, and returns a
polynomial representing the `i,j` th entry of the ``round`` th round key. 
This polynomial's variables are entries of the various subkeys from a key
expansion. Note, however, that while all polynomials with key variables will
have variables from subkey entries, the variables of the polynomials returned
by ``expand_key_poly`` will be in terms of

    sage: rgf.expand_key_poly(1, 2, 0)
    k012
    sage: rgf.expand_key_poly(1, 1, 1)
    k111
    sage: rgf.expand_key_poly(1, 2, 1)
    (x^2 + 1)*k121^254 + (x^3 + 1)*k121^253 + (x^7 + x^6 + x^5 + x^4 + x^3 + 1)*k121^251 + (x^5 + x^2 + 1)*k121^247 + (x^7 + x^6 + x^5 + x^4 + x^2)*k121^239 + k121^223 + (x^7 + x^5 + x^4 + x^2 + 1)*k121^191 + (x^7 + x^3 + x^2 + x + 1)*k121^127 + k010 + (x^6 + x^5 + x)

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
from sage.structure.element import Matrix
from sage.rings.finite_rings.constructor import FiniteField
from sage.rings.integer import Integer
from sage.structure.sage_object import SageObject
from sage.matrix.matrix_space import MatrixSpace
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing

class RijndaelGF(SageObject):

    def __init__(self, Nb, Nk, state_chr='a', key_chr='k'):
        r"""
        An algebraically generalized version of the AES cipher.
        
        INPUT:
        
        - ``Nb`` -- The block length of this instantiation. Must be between 4
          and 8.
        
        - ``Nk`` -- The key length of this instantion. Must be between 4 and 8.

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

        We can also alter the name of variables in polynomials representing
        elements from round keys by changing ``key_chr``. ::

            sage: rgf = RijndaelGF(4, 6, key_chr='myKeyChr')
            sage: rgf.expand_key_poly(1, 2, 1)
            (x^2 + 1)*myKeyChr121^254 + (x^3 + 1)*myKeyChr121^253 + (x^7 + x^6 + x^5 + x^4 + x^3 + 1)*myKeyChr121^251 + (x^5 + x^2 + 1)*myKeyChr121^247 + (x^7 + x^6 + x^5 + x^4 + x^2)*myKeyChr121^239 + myKeyChr121^223 + (x^7 + x^5 + x^4 + x^2 + 1)*myKeyChr121^191 + (x^7 + x^3 + x^2 + x + 1)*myKeyChr121^127 + myKeyChr010 + (x^6 + x^5 + x)
        """
        if Nb not in range(4, 9):
            msg = "Block length Nb must be in the range 4 - 8, not {0}"
            raise ValueError(msg.format(Nb))
        if Nk not in range(4, 9):
            msg = "Key length Nk must be in the range 4 - 8, not {0}"
            raise ValueError(msg.format(Nk))
        if not isinstance(state_chr, basestring):
            msg = "state_chr must be a string, not {0}"
            raise TypeError(msg.format(state_chr))
        if not isinstance(key_chr, basestring):
            msg = "key_chr must be a string, not {0}"
            raise TypeError(msg.format(key_chr))

        self._Nb = Nb
        self._Nk = Nk
        from sage.rings.polynomial.polynomial_ring import polygen
        from sage.rings.finite_rings.integer_mod_ring import Integers
        pgen = polygen(Integers(2))
        mod = pgen**8 + pgen**4 + pgen**3 + pgen + 1
        self._F = FiniteField(2**8, 'x', modulus=mod)
        round_num_table = matrix([[10,11,12,13,14], [11,11,12,13,14],
                                  [12,12,12,13,14], [13,13,13,13,14],
                                  [14,14,14,14,14]])
        self._Nr = round_num_table[self._Nb - 4, self._Nk - 4]

        # Build framework for polynomial creation.
        state_names = [state_chr + str(i) + str(j) 
                       for i in range(4) for j in range(self._Nb)]
        subkey_names = [key_chr + str(r) + str(i) + str(j)
                        for r in range(self._Nr + 1) for i in range(4)
                        for j in range(self._Nb)]
        self._state_PR = PolynomialRing(self._F, len(state_names), state_names)
        self._all_PR = PolynomialRing(self._F, len(state_names + subkey_names),
                                      state_names + subkey_names)
        self.state_vrs = matrix(4, self._Nb, self._state_PR.gens())
        self.subkey_vrs_list = list(self._all_PR.gens()[4 * self._Nb:])
        self.subkey_vrs = [matrix(4, self._Nb,
                           self.subkey_vrs_list[(4 * self._Nb)*i :
                                        (4 * self._Nb)*(i+1)])
                           for i in range(self._Nr)]
        self.key_vrs = column_matrix([
                       self.subkey_vrs[int(i / self._Nb)].column(i % 4)
                       for i in range(self._Nk)])
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
        Returns the encryption/decryption of ``text`` with key ``key``.

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

        We can also use binary strings for ``text`` and ``key``. ::

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
        Returns the string representation of ``self``.

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

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF 
            sage: rgf = RijndaelGF(5, 4)
            sage: rgf.number_rounds()
            11
        """
        return self._Nr
        
    def hex_to_GF(self, H, matrix=True):
        r"""
        Returns a matrix/list of elements of `\GF{2^8}` corresponding to ``H``.

        INPUT:
        
        - ``H`` -- A hex string where every two hex characters correspond to a
          single element in `\GF{2^8}`
          
        - ``matrix`` -- (default: ``True``) Returns a list if ``False``;
          returns a state matrix if ``True``.

        OUTPUT:

        - A list of or a state matrix of elements of `\GF{2^8}` where each
          element corresponds to the appropriate hex value in ``H``. In
          particular, every element `a_7x^7 + a_6x^6 + a_5x^5 + a_4x^4 +
          a_3x^3 + a_2x^2 + a_1x^1 + a_0` in `\GF{2^8}` corresponds to the
          8-bit binary string '`a_7a_6a_5a_4a_3a_2a_1a_0`'.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: state = rgf.hex_to_GF('1147659047cf663b9b0ece8dfc0bf1f0')
            sage: output = rgf.shift_rows(state)
            sage: rgf.GF_to_hex(output)
            '11cfcef0470ef1909b0b653bfc47668d'

        We can output a list instead by setting ``matrix`` to ``False``. ::

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
        Returns the hex string representation of ``GF``.

        INPUT:
        
        - ``GF`` -- Either a state matrix over `\GF{2^8}`, a list of elements
          from `\GF{2^8}`, or a single element from `\GF{2^8}`

        OUTPUT:

        - A hex string representation of ``GF``, where every two characters in
          the string correspond to a single element in `\GF{2^8}`. In
          particular, every element `a_7x^7 + a_6x^6 + a_5x^5 + a_4x^4 +
          a_3x^3 + a_2x^2 + a_1x^1 + a_0` in `\GF{2^8}` corresponds to the
          8-bit binary string '`a_7a_6a_5a_4a_3a_2a_1a_0`'.

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
        if not isinstance(GF, Matrix) and \
           not isinstance(GF, list) and \
           not is_FiniteFieldElement(GF):
            msg = ("keyword 'GF' must be a matrix over {0}, a list of "
                   "elements from {0}, or a single element from {0}")
            raise TypeError(msg.format(self._F))

        if isinstance(GF, Matrix):
            if not GF.base_ring().is_field() or \
               not GF.base_ring().is_finite() or \
               not GF.base_ring().order() == 2**8:
                msg = "The elements of keyword 'GF' must all be from {0}"
                raise TypeError(msg.format(self._F))
            return ''.join([self.GF_to_hex(el)
                            for col in GF.columns() for el in col])
        elif isinstance(GF, list):
            if not all([g.parent().is_field() and g.parent().is_finite() and
                        g.parent().order() == 2**8 for g in GF]):
                msg = "The elements of keyword 'GF' must all be from {0}"
                raise TypeError(msg.format(self._F))
            return ''.join([self.GF_to_hex(el) for el in GF])
        else:
            if not GF.parent().is_field() or \
               not GF.parent().is_finite() or \
               not GF.parent().order() == 2**8:
                msg = "keyword 'GF' must be in"
                raise TypeError(msg.format(self._F))
            return hex(GF.integer_representation())[2:].zfill(2)

    def bin_to_GF(self, B, matrix=True):
        r"""
        Returns a matrix/list of elements of `\GF{2^8}` corresponding to ``B``.

        INPUT:

        - ``B`` -- A binary string where every eight bits correspond to a
          single element in `\GF{2^8}`
          
        - ``matrix`` -- (default: ``True``) Returns a list if ``False``.
          Returns a state matrix over `\GF{2^8}` if ``True``.

        OUTPUT:

        - A list of or a state matrix of elements of `\GF{2^8}` where each
          element corresponds to the appropriate 8-bit binary string in ``B``.
          In particular, every element `a_7x^7 + a_6x^6 + a_5x^5 + a_4x^4 +
          a_3x^3 + a_2x^2 + a_1x^1 + a_0` in `\GF{2^8}` corresponds to the
          8-bit binary string '`a_7a_6a_5a_4a_3a_2a_1a_0`'.

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
        Returns the binary string representation of ``GF``.

        INPUT:
        
        - ``GF`` -- Either a state matrix over `\GF{2^8}`, a list of elements
          from `\GF{2^8}`, or a single element from `\GF{2^8}`

        OUTPUT:

        - A binary string representation of ``GF``, where every eight 
          characters in the string corresponds to a single element in
          `\GF{2^8}`. In particular, every element `a_7x^7 + a_6x^6 + a_5x^5 +
          a_4x^4 + a_3x^3 + a_2x^2 + a_1x^1 + a_0` in `\GF{2^8}` corresponds
          to the binary string '`a_7a_6a_5a_4a_3a_2a_1a_0`'.

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
        if not isinstance(GF, Matrix) and \
           not isinstance(GF, list) and \
           not is_FiniteFieldElement(GF):
            msg = ("keyword 'GF' must be a matrix over {0}, a list of "
                   "elements from {0}, or a single element from {0}")
            raise TypeError(msg.format(self))

        if isinstance(GF, Matrix):
            if not GF.base_ring().is_field() or \
               not GF.base_ring().is_finite() or \
               not GF.base_ring().order() == 2**8:
                msg = "The elements of keyword 'GF' must all be from {0}"
                raise TypeError(msg.format(self._F))
            return ''.join([self.GF_to_bin(el)
                            for col in GF.columns() for el in col])
        elif isinstance(GF, list):
            if not all([g.parent().is_field() and g.parent().is_finite() and
                        g.parent().order() == 2**8 for g in GF]):
                msg = "The elements of keyword 'GF' must all be from {0}"
                raise TypeError(msg.format(self._F))
            return ''.join([self.GF_to_bin(el) for el in GF])
        else:
            if not GF.parent().is_field() or \
               not GF.parent().is_finite() or \
               not GF.parent().order() == 2**8:
                msg = "keyword 'GF' must be in"
                raise TypeError(msg.format(self._F))
            return bin(GF.integer_representation())[2:].zfill(8)

    def encrypt(self, plain, key, format='hex'):
        r"""
        Returns the plaintext ``plain`` encrypted with the key ``key``.

        INPUT:

        - ``plain`` -- The plaintext to be encrypted.

        - ``key`` -- The key to encrypt ``plain`` with.

        - ``format`` -- (default: ``hex``) The string format of ``key`` and
          ``plain``, either "hex" or "binary".

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
        Returns the ciphertext ``ciphertext`` decrypted with the key ``key``.

        INPUT:

        - ``ciphertext`` -- The ciphertext to be decrypted.
        
        - ``key`` -- The key to decrypt ``ciphertext`` with.

        - ``format`` -- (default: ``hex``) The string format that both
          ``ciphertext`` and ``key`` must be in, either "hex" or "binary".

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
                  
    def _check_valid_PRmatrix(self, PRm, keyword):
        r"""
        Raises an error if ``PRm`` is not a valid input matrix.

        INPUT:

        - ``PRm`` -- If ``PRm`` is a `4 \times Nb` matrix with entries from
          the multivariate PolynomialRing ``_all_PR``, this method does nothing
          `\GF{2^8}`, this method does nothing. Otherwise, this method raises
          an error. Note that a matrix of elements from `\GF(2^8)` is regarded
          as a matrix with entries from ``_all_PR`` and will pass this test.

        - ``keyword`` -- The name of the keyword ``PRm`` from where this
          method was called, for the potential error message. For example, if
          called from ``sub_bytes``, ``keyword`` would be "state".
          
        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: good_state = rgf.hex_to_GF('0'*32)
            sage: rgf._check_valid_PRmatrix(good_state, 'state')
            sage: rgf._check_valid_PRmatrix(rgf.state_vrs, 'state')
            sage: rgf._check_valid_PRmatrix(5, 'state')
            Traceback (most recent call last):
            ...
            TypeError: keyword 'state' must be a 4 x 4 matrix with entries from a multivariate PolynomialRing over Finite Field in x of size 2^8
            <BLANKLINE>
            sage: entries = [rgf._F.random_element() for i in range(24)]
            sage: wrong_dimensions = matrix(4, 6, entries)
            sage: rgf._check_valid_PRmatrix(wrong_dimensions, 'state')
            Traceback (most recent call last):
            ...
            TypeError: keyword 'state' must be a 4 x 4 matrix with entries from a multivariate PolynomialRing over Finite Field in x of size 2^8
            <BLANKLINE>
            sage: F.<a> = GF(3^4)
            sage: entries = [F.random_element() for i in range(16)]
            sage: wrong_base = matrix(4, 4, entries)
            sage: rgf._check_valid_PRmatrix(wrong_base, 'state')
            Traceback (most recent call last):
            ...
            TypeError: keyword 'state' must be a 4 x 4 matrix with entries from a multivariate PolynomialRing over Finite Field in x of size 2^8
        """
        from sage.rings.polynomial.multi_polynomial_ring_generic import \
            MPolynomialRing_generic
        msg = ("keyword '{0}' must be a {1} x {2} matrix with entries from a "
               "multivariate PolynomialRing over {3}")
        msg = msg.format(keyword, 4, self._Nb, self._F)
        if (not isinstance(PRm, Matrix) or \
            not (PRm.base_ring().is_field() and \
                PRm.base_ring().is_finite() and \
                PRm.base_ring().order() == 256 and \
                PRm.dimensions() == (4, self._Nb))) and \
           (not isinstance(PRm, Matrix) or \
            not isinstance(PRm.base_ring(), MPolynomialRing_generic) or \
            not (PRm.base_ring().base_ring().is_field() and \
                 PRm.base_ring().base_ring().is_finite() and \
                 PRm.base_ring().base_ring().order() == 256) or \
                not PRm.dimensions() == (4, self._Nb)):
            raise TypeError(msg)

    def _test_poly_input_bounds(self, row, col):
        r"""
        Raises an error if ``[row,col]`` is an invalid state matrix entry.

        INPUT: 

        - ``row`` -- The row index.

        - ``col`` -- The column index.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(6, 4)            
            sage: rgf._test_poly_input_bounds(0, 0)
            sage: rgf._test_poly_input_bounds(5, 0)
            Traceback (most recent call last): 
            ...
            ValueError: keyword 'row' must be between 0 and 4
            sage: rgf._test_poly_input_bounds(0, 7)
            Traceback (most recent call last):
            ...
            ValueError: keyword 'col' must be between 0 and 6
        """
        if row not in range(0, 4):
            raise ValueError("keyword 'row' must be between 0 and 4")
        if col not in range(0, self._Nb):
            msg = "keyword 'col' must be between 0 and {0}"
            raise ValueError(msg.format(self._Nb))

    def expand_key(self, key):
        r"""
        Returns the expanded key schedule from ``key``.

        INPUT:

        - ``key`` -- The key to build a key schedule from. Must be a matrix
          over `\GF{2^8}` of dimensions `4 \times N_k`.

        OUTPUT:

        - A length `Nr` list of `4 \times N_b` matrices corresponding to the
          expanded key. The `n` th entry of the list corresponds to the matrix
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
        msg = "keyword '{0}' must be a {1} x {2} matrix over GF({3})"
        msg = msg.format(key, 4, self._Nk, self._F.order())
        if not isinstance(key, Matrix) or \
           not (key.base_ring().is_field() and \
                key.base_ring().is_finite() and \
                key.base_ring().order() == self._F.order()) or \
           not key.dimensions() == (4, self._Nk):
            raise TypeError(msg)

        def add_cols(col1, col2):
            return map(lambda (x,y): x + y, zip(col1, col2))

        key_cols = []
        for i in range(self._Nb * (self._Nr + 1)):
            key_cols.append([])

        # Copy columns from ``key``, then build the rest of the columns
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

        # Copy the expanded columns into 4xNb blocks
        round_keys = []
        for r in range(self._Nr + 1):
            rk = column_matrix([key_cols[r*self._Nb + i] 
                                for i in range(self._Nb)])
            round_keys.append(rk)
        return round_keys

    def expand_key_poly(self, row, col, round):
        r"""
        Returns a polynomial representing the ``row,col`` th entry of the
        ``round`` th round key.

        INPUT:

        - ``row`` -- The row position of the element represented by this
          polynomial.
          
        - ``col`` -- The column position of the element represented by this
          polynomial.

        OUTPUT:

        - A polynomial representing the ``row,col`` th entry of the ``round``
          th round key in terms of entries of the input key.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: rgf.expand_key_poly(1, 2, 0)
            k012
            sage: rgf.expand_key_poly(1, 2, 1)
            (x^2 + 1)*k023^254 + (x^3 + 1)*k023^253 + (x^7 + x^6 + x^5 + x^4 + x^3 + 1)*k023^251 + (x^5 + x^2 + 1)*k023^247 + (x^7 + x^6 + x^5 + x^4 + x^2)*k023^239 + k023^223 + (x^7 + x^5 + x^4 + x^2 + 1)*k023^191 + (x^7 + x^3 + x^2 + x + 1)*k023^127 + k010 + k011 + k012 + (x^6 + x^5 + x)

        It should be noted that ``expand_key_poly`` cannot be used with
        ``apply_poly`` or ``compose``, due to the state polynomials and the
        expanded key polynomials being represented differently. ::

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
        self._test_poly_input_bounds(row, col)
        if round not in range(self._Nr + 1):
            msg = "keyword 'r' must be between 0 and {0}"
            raise ValueError(msg.format(self._Nr))

        key_col = round * self._Nb + col
        if key_col < self._Nk:
            return self.key_vrs[row, key_col]
        else:

            if key_col % self._Nk == 0 or \
               (self._Nk > 6 and col % self._Nk == 4):
                # Apply non-linear transformation to key_col - 1
                recur_r = int((key_col - 1)/self._Nb)
                recur_j = (key_col - 1) - (recur_r * self._Nb)
                non_linear = self.expand_key_poly((row+1) % 4, 
                                                  recur_j, recur_r)
                non_linear = self._srd(non_linear)
                non_linear += self._F.gen() ** (int(key_col / self._Nk) - 1)
                # Identify key_col - Nk
                recur_r = int((key_col - self._Nk)/self._Nb)
                recur_j = (key_col - self._Nk) - (recur_r * self._Nb)
                return self.expand_key_poly(row, recur_j, recur_r) + non_linear
            else:
                # Identify key_col - Nk
                recur_r = int((key_col - self._Nk)/self._Nb)
                recur_j = (key_col - self._Nk) - (recur_r * self._Nb)
                result = self.expand_key_poly(row, recur_j, recur_r)
                # Identify key_col - 1
                recur_r = int((key_col- 1)/self._Nb)
                recur_j = (key_col - 1) - (recur_r * self._Nb)
                return result + self.expand_key_poly(row, recur_j, recur_r)
        
    def apply_poly(self, state, poly_method, algorithm='encrypt', keys=None):
        r"""
        Returns a state matrix where ``poly_method`` is applied to each entry.

        INPUT:
        
        - ``state`` -- The state matrix over `\GF{2^8}` to which 
          ``poly_method`` is applied to.

        - ``poly_method`` -- The method to apply to ``state``. ``poly_method``
          must be callable as ``poly_method(i, j, algorithm)`` and must return
          a polynomial.

        - ``algorithm`` -- (default: "encrypt") Passed directly to 
          ``poly_method`` to select encryption or decryption. The
          encryption flag is "encrypt" and the decrypt flag is "decrypt".

        - ``keys`` -- (default: None) An array of `N_r` subkey matrices to
          replace any key variables in any polynomials returned by
          ``poly_method``. Must be identical to the format returned by
          ``expand_key``. If any polynomials have key variables and ``keys``
          is not supplied, the key variables will remain as-is.

        OUTPUT:

        - A state matrix in `\GF{2^8}` whose `i,j` th entry equals the
          polynomial ``poly_method(i, j, algorithm)`` evaluated by setting
          its variables equal to the corresponding entries of ``state``.

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

        If ``poly_method`` returns a polynomial with state variables as well
        as key variables, we can supply a list of round keys ``keys`` whose
        elements are evaluated as the key variables. ::

            sage: def fn(row, col, algorithm='encrypt'):
            ....:     return rgf.add_round_key_poly(row, col, 0, algorithm)
            ....:
            sage: fn(1, 2)
            a12 + k012
            sage: state = rgf.hex_to_GF('14f9701ae35fe28c440adf4d4ea9c026')
            sage: key = rgf.hex_to_GF('54d990a16ba09ab596bbf40ea111702f')
            sage: keys = rgf.expand_key(key)
            sage: result = rgf.apply_poly(state, fn, keys)
        """
        self._check_valid_PRmatrix(state, 'state')
        if poly_method == self.expand_key_poly:
            raise ValueError("expand_key_poly cannot be used with apply_poly")
        if keys != None and (not isinstance(keys, list) or \
           len(keys) != self._Nr + 1 or \
           not all([isinstance(k, Matrix) for k in keys]) or \
           not all([k.dimensions() == (4, self._Nb) for k in keys]) or \
           not all([k.base_ring().is_finite() and k.base_ring().is_field()
                    and k.base_ring().order() == 256 for k in keys]) ):
            msg = ("keys must be a length {0} array of 4 by {1} matrices"
                   " over {2}")
            raise TypeError(msg.format(self._Nr, self._Nb, self._F))

        output = []
        if keys != None:
            key_list = [el for inner in keys for el in inner.list()]
        for i in range(4):
            for j in range(self._Nb):
                # this is to combat a major performance issue caused by 
                # subbytes' inversion transformation.
                if poly_method == self.sub_bytes_poly and \
                   algorithm == 'decrypt':
                    p = poly_method(i, j, algorithm, no_inversion=True)
                    p = p(state.list()) ** 254
                else:
                    p = poly_method(i, j, algorithm)
                    # If there are key variables in the polynomial
                    if len(p.args()) > 4 * self._Nb:
                        if keys != None:
                            p = p(state.list() + key_list)
                        else:
                            p = p(state.list() + self.subkey_vrs_list)
                    else:
                        p = p(state.list())
                output.append(p)
        return matrix(4, 4, output)

    def compose(self, f, g, algorithm='encrypt'):
        r"""
        Returns a polynomial representing `g(f(A))_{i,j}` if ``g`` is a
        polynomial representing `g(A)_{i,j}`, or returns a ``_poly`` method
        which returns polynomials of that form if `g` is a ``_poly`` method.

        INPUT:

        - ``f`` -- A ``_poly`` function of the form ``f(i, j, algorithm)``
          which returns a polynomial reprentation of `f(A)_{i,j}`.

        - ``g`` -- A ``_poly`` function of the form ``g(i, j, algorithm)``
          which returns a polynomial reprentation of `g(A)_{i,j}` OR a
          polynomial representing `g(A)_{i,j}`.

        - ``algorithm`` -- (default: "encrypt") Whether ``f`` and ``g``
          should use their encryption transformations or their decryption
          transformations. Does nothing if ``g`` is a function. The encryption
          flag is "encrypt" and the decryption flag is "decrypt".

        OUTPUT:

        - If ``g`` is a method which takes an index ``i,j`` and returns a
          polynomial representing the `i,j` th entry of the matrix `g(A)` in
          terms of entries of `A`, where `A` is a state matrix, then
          ``compose(f,g)`` returns a method which takes an index ``i,j`` and
          returns a polynomial representing the `i,j` th entry of the matrix
          `g(f(A))` in terms of entries of `A`. On the other hand, if ``g`` is
          a polynomial representing the `i,j` th entry of the matrix `g(A)` in
          terms of entries of `A`, then ``compose(f,g)`` returns a polynomial
          representing the `i,j` th entry of `g(f(A))` in terms of entries of
          `A`. In both circumstances, ``f`` must be a method which takes an
          index ``i,j`` and returns a polynomial representing the `i,j` th
          entry of `f(A)` in terms of the entries of the state matrix `A`.

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

        If the second argument is a polynomial, then the value of ``algorithm``
        is passed directly to the first argument `f` during evaluation.
        However, if the second argument is a method, changing ``algorithm``
        does nothing since the returned method must have its own settable
        ``algorithm`` keyword in order to have the same arguments as the
        rest of the ``_poly(row, col, algorithm='encrypt')`` methods. ::
        
            sage: f = rgf.compose(rgf.sub_bytes_poly, rgf.mix_columns_poly,
            ....: algorithm='decrypt')
            sage: g = rgf.compose(rgf.sub_bytes_poly, rgf.mix_columns_poly)
            sage: all([f(i,j) == g(i,j) for i in range(4) for j in range(4)])
            True
        """
        if f == self.expand_key_poly or g == self.expand_key_poly:
            raise ValueError("expand_key_poly cannot be used with compose")
        if g in self._all_PR:
            f_vals = [f(i, j, algorithm)
                      for i in range(4) for j in range(self._Nb)]
            if g in self._state_PR:
                return g(f_vals)
            else:
                return g(f_vals + self.subkey_vrs_list)
        else:
            lm = lambda i, j, alg='encrypt': self.compose(f, g(i, j, alg), alg)
            return lm

    def add_round_key_poly(self, row, col, round, algorithm='encrypt'):
        r"""
        Returns a polynomial representing an element of a round-key addition.

        INPUT:

        - ``row`` -- The row number of the entry represented by this method's
          output.

        - ``col`` -- The column number of the entry represented by this
          method's output.

        - ``round`` -- The round number of the entry represented by this
          method's output.

        - ``algorithm`` -- (default: "encrypt") Whether to return the
          polynomial as an encryption or as a decryption. The encryption flag
          is "encrypt" and the decryption  flag is "decrypt".

        OUTPUT:

        - A polynomial representing the ``row,col`` th entry of a state matrix
          after a round-key addition in terms of entries of the input state
          matrix and entries of the ``round`` th round key.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: rgf.add_round_key_poly(1, 2, 7)
            a12 + k712

        As expected, since the encryption and decryption transformations are
        identical, changing ``algorithm`` has no effect.

            sage: with_encrypt = rgf.add_round_key_poly(3, 2, 6, 'encrypt')
            sage: with_decrypt = rgf.add_round_key_poly(3, 2, 6, 'decrypt')
            sage: with_encrypt == with_decrypt
            True
        """
        self._test_poly_input_bounds(row, col)
        if round not in range(self._Nr):
            msg = "keyword 'round' must be between 0 and {0}"
            raise ValueError(msg.format(self._Nr))
        state_var = self.state_vrs[row, col]
        key_var = self.subkey_vrs[round][row, col]
        return state_var + key_var

    def add_round_key(self, state, round_key):
        r"""
        Returns the round-key addition of matrices ``state`` and ``round_key``.

        INPUT:
        
        - ``state`` -- The state matrix to have ``round_key`` added to.

        - ``round_key`` -- The round key to add to ``state``.

        OUTPUT:

        - A state matrix which is the round key addition of ``state`` and
          ``round_key``. This transformation is simply the entrywise addition
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
        self._check_valid_PRmatrix(state, 'state')
        self._check_valid_PRmatrix(round_key, 'round_key')
        # We don't use apply_poly here since that would require an extra
        # argument of the round number
        return state + round_key

    def sub_bytes_poly(self, row, col, algorithm='encrypt', 
                       no_inversion=False):
        r"""
        Returns a polynomial representing `SubBytes(A)_{\textit{row, col}}`.

        INPUT:

        - ``row`` -- The row number of the entry represented by this method's
          output.

        - ``col`` -- The column number of the entry represented by this
          method's output.

        - ``algorithm`` -- (default: "encrypt") Whether to return the
          polynomial as an encryption or as a decryption. The encryption flag
          is "encrypt" and the decryption  flag is "decrypt".

        - ``no_inversion`` -- (default: ``False``) Don't perform the inversion
          step, only perform the affine transformation. Primarily intended
          to increase performance during decryption, as is shown in the
          below example.

        OUTPUT:

        - A polynomial representing the ``row,col`` th entry of a state matrix
          after the SubBytes method has been applied to it.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: rgf.sub_bytes_poly(2, 3)
            (x^2 + 1)*a23^254 + (x^3 + 1)*a23^253 + (x^7 + x^6 + x^5 + x^4 + x^3 + 1)*a23^251 + (x^5 + x^2 + 1)*a23^247 + (x^7 + x^6 + x^5 + x^4 + x^2)*a23^239 + a23^223 + (x^7 + x^5 + x^4 + x^2 + 1)*a23^191 + (x^7 + x^3 + x^2 + x + 1)*a23^127 + (x^6 + x^5 + x + 1)

        We can use this polynomial to calculate individual entries of the
        output matrix for any given state as such: ::

            sage: state = rgf.hex_to_GF('6385b79ffc538df997be478e7547d691')
            sage: poly = rgf.sub_bytes_poly(2, 3)
            sage: poly(state.list())
            x^7 + x^6 + x^5 + x^4 + x^2 + x

        We can set ``no_inversion`` to ``True`` to get a polynomial
        representation of solely the affine transformation. ::

            sage: rgf.sub_bytes_poly(0, 2, no_inversion=True)
            (x^7 + x^3 + x^2 + x + 1)*a02^128 + (x^7 + x^5 + x^4 + x^2 + 1)*a02^64 + a02^32 + (x^7 + x^6 + x^5 + x^4 + x^2)*a02^16 + (x^5 + x^2 + 1)*a02^8 + (x^7 + x^6 + x^5 + x^4 + x^3 + 1)*a02^4 + (x^3 + 1)*a02^2 + (x^2 + 1)*a02 + (x^6 + x^5 + x + 1)

        When generating a decryption polynomial, calculating the inverse of
        the polynomial representing the affine transformation can be a very
        slow process. In order to speed up decryption when applying
        ``sub_bytes_poly`` to a state matrix, it is recommended to calculate
        the decryption polynomial with ``no_inversion=True``, evaluate the
        arguments, then perform the inversion after this result has been
        calculated. ::

            sage: poly = rgf.sub_bytes_poly(0, 0, 
            ....: algorithm='decrypt', no_inversion=True)
            sage: state = rgf.hex_to_GF('b415f8016858552e4bb6124c5f998a4c')
            sage: poly(state.list()) ^ -1
            x^7 + x^6 + x^2 + x
        """
        self._test_poly_input_bounds(row, col)
        if algorithm == 'encrypt':
            var = self.state_vrs[row, col]
            coeffs = self._sb_E_coeffs
            if no_inversion:
                return sum([coeffs[i] * (var**(2**i)) 
                            for i in range(8)]) + self._F("x^6 + x^5 + x + 1")
            else:
                return sum([coeffs[i] * (var**(255 - 2**i)) 
                            for i in range(8)]) + self._F("x^6 + x^5 + x + 1")
        elif algorithm == 'decrypt':
            var = self.state_vrs[row, col]
            coeffs = self._sb_D_coeffs
            result = (sum([coeffs[i] * var**(2**i) for i in range(8)]) + \
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
        Returns the application of SubBytes (`S_{RD}`) to ``el``.

        INPUT:

        - ``el`` -- An element of `\GF{2^8}`.

        - ``algorithm`` -- (default: "encrypt") Whether to perform the
          encryption transformation or the decryption transformation.
          The encryption flag is "encrypt" and the decryption flag is
          "decrypt".

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
            state = [el] + [self._F.zero()]*((4 * self._Nb)-1)
            return p(state)
        elif algorithm == 'decrypt':
            p = self.sub_bytes_poly(0, 0, algorithm, no_inversion=True)
            state = [el] + [self._F.zero()]*((4 * self._Nb)-1)
            return p(state) ** 254
        else:
            raise ValueError(("keyword 'algorithm' must be either 'encrypt' "
                             "or 'decrypt'"))

    def sub_bytes(self, state, algorithm='encrypt'):
        r"""
        Returns the application of `SubBytes` to the state matrix ``state``.

        INPUT:

        - ``state`` -- The state matrix to apply SubBytes to.
        
        - ``algorithm`` -- (default: "encrypt") Whether to apply the
          encryption step of SubBytes or its decryption inverse. The encryption
          flag is "encrypt" and the decryption flag is "decrypt".

        OUTPUT:

        - The state matrix over `\GF{2^8}` where `SubBytes` has been applied
          to every entry of ``state``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: state = rgf.hex_to_GF('d1c4941f7955f40fb46f6c0ad68730ad')
            sage: result = rgf.sub_bytes(state)
            sage: rgf.GF_to_hex(result)
            '3e1c22c0b6fcbf768da85067f6170495'
            sage: decryption = rgf.sub_bytes(result, algorithm='decrypt')
            sage: decryption == state
            True
        """
        self._check_valid_PRmatrix(state, 'state')
        return self.apply_poly(state, self.sub_bytes_poly, algorithm)

    def mix_columns_poly(self, row, col, algorithm='encrypt'):
        r"""
        Returns a polynomial representing `MixColumns(A)_{\textit{row, col}}`.

        INPUT:
        
        - ``row`` -- The row number of the entry represented by this method's
          output.
          
        - ``col`` -- The column number of the entry represented by this
          method's output.

        - ``algorithm`` -- (default: "encrypt") Whether to perform the
          encryption transformation or the decryption transformation. The
          encryption flag is "encrypt" and the decryption flag is "decrypt".

        OUTPUT:
        
        - A polynomial in terms of entries of the input state matrix which
          represents the ``row,col`` th entry of the output matrix after
          MixColumns has been applied to it.

        EXAMPLES::
        
            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: rgf.mix_columns_poly(3, 1)
            (x + 1)*a01 + a11 + a21 + (x)*a31

        We can use this to calculate individual entries of a state matrix after
        the decryption version of MixColumns has been applied to it as such: ::

            sage: poly = rgf.mix_columns_poly(2, 2, algorithm='decrypt')
            sage: state = rgf.hex_to_GF('a761ca9b97be8b45d8ad1a611fc97369')
            sage: result = poly(state.list())
            sage: rgf.GF_to_hex(result)
            'b7'
            sage: output = rgf.mix_columns(state, algorithm='decrypt')
            sage: output[2,2] == result
            True
        """
        self._test_poly_input_bounds(row, col)
        if algorithm == 'encrypt':
            coeffs = self._mixcols_E
        elif algorithm == 'decrypt':
            coeffs = self._mixcols_D
        else:
            raise ValueError(("keyword 'algorithm' must be either 'encrypt' "
                             "or 'decrypt'"))
        return sum([coeffs[row,k] * self.state_vrs[k,col] for k in range(4)])

    def mix_columns(self, state, algorithm='encrypt'):
        r"""
        Returns the application of `MixColumns` to the state matrix ``state``.

        INPUT:

        - ``state`` -- The state matrix to apply MixColumns to.
        
        - ``algorithm`` -- (default: "encrypt") Whether to perform the
          encryption version of MixColumns, or its decryption inverse. The
          encryption flag is "encrypt" and the decryption flag is "decrypt".

        OUTPUT:
        
        - The state matrix over `\GF{2^8}` which is the result of applying
          MixColumns to ``state``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: state = rgf.hex_to_GF('cd54c7283864c0c55d4c727e90c9a465')
            sage: result = rgf.mix_columns(state)
            sage: rgf.GF_to_hex(result)
            '921f748fd96e937d622d7725ba8ba50c'
            sage: decryption = rgf.mix_columns(result, algorithm='decrypt')
            sage: decryption == state
            True
        """
        self._check_valid_PRmatrix(state, 'state')
        return self.apply_poly(state, self.mix_columns_poly, algorithm)

    def shift_rows_poly(self, row, col, algorithm='encrypt'):
        r"""
        Returns a polynomial representing `ShiftRows(A)_{\textit{row,col}}`.

        INPUT:

        - ``row`` -- The row number of the entry represented by this method's
          output.

        - ``col`` -- The column number of the entry represented by this
          method's output.

        - ``algorithm`` -- (default: "encrypt") Whether to perform ShiftRows'
          encryption step or its decryption inverse. The encryption flag is
          "encrypt" and the decryption flag is "decrypt".

        OUTPUT:
        
        - A polynomial in terms of entries of the input state matrix which
          represents the ``row,col`` th entry of the output matrix after
          ShiftRows has been applied to it.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4, 4)
            sage: rgf.shift_rows_poly(2, 3)
            a21
            
        We can use this to calculate individual entries of a state matrix
        after the decryption version of ShiftRows has been applied to it as
        such: ::

            sage: poly = rgf.shift_rows_poly(2, 3, algorithm='decrypt')
            sage: state = rgf.hex_to_GF('78c4f708318d3cd69655b701bfc093cf')
            sage: result = poly(state.list())
            sage: rgf.GF_to_hex(result)
            '3c'
            sage: output = rgf.shift_rows(state, algorithm='decrypt')
            sage: output[2,3] == result
            True        
        """
        self._test_poly_input_bounds(row, col)
        if algorithm == 'encrypt':
            offs = self._shiftrows_offsets_E
        elif algorithm == 'decrypt':
            offs = self._shiftrows_offsets_D
        else:
            raise ValueError(("keyword 'algorithm' must be either 'encrypt' "
                             "or 'decrypt'"))
        return self.state_vrs[row, (col + offs[4 - self._Nb][row]) % 4]
    
    def shift_rows(self, state, algorithm='encrypt'):
        r"""
        Returns the application of `ShiftRows` to the state matrix ``state``.

        INPUT:
        
        - ``state`` -- A state matrix over `\GF{2^8}` to which ShiftRows is 
          applied to.

        - ``algorithm`` -- (default: "encrypt") Whether to perform the 
          encryption version of ShiftRows or its decryption inverse. The 
          encryption flag is "encrypt" and the decryption flag is "decrypt".

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
            sage: decryption = rgf.shift_rows(result, algorithm='decrypt')
            sage: decryption == state
            True
        """
        self._check_valid_PRmatrix(state, 'state')
        return self.apply_poly(state, self.shift_rows_poly, algorithm)
