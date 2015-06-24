r"""
Rijndael-GF

An algebraic implementation of the Rijndael Cipher (more commonly 
known as AES). This implementation is generalized to allow for an algebraic
treatment of the cipher and the individual steps of the round for the purposes
of learning, study in relation to other algebraic ciphers, and algebraic
cryptanalysis.

AUTHORS:

- Thomas Gagne (2015-06): initial version
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
    This class implements the Rijndael-GF cipher as described in [RIJNDAEL].
    """
    
    def __init__(self, Nb, Nk):
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
        self._shiftrows_offsets_E = matrix([[0,1,2,3], [0,1,2,3], [0,1,2,3],
                                   [0,1,2,4], [0,1,3,4]])
        self._shiftrows_offsets_D = matrix([[0,-1,-2,-3], [0,-1,-2,-3],
                                           [0,-1,-2,-3], [0,-1,-2,-4],
                                           [0,-1,-3,-4]])
        self._polyring = PolynomialRing(self._F, 'p')
        sb_E_coeffs = [self._F("x^2 + 1"), 
                       self._F("x^3 + 1"), 
                       self._F("x^7 + x^6 + x^5 + x^4 + x^3 + 1"), 
                       self._F("x^5 + x^2 + 1"), 
                       self._F("x^7 + x^6 + x^5 + x^4 + x^2"), 
                       self._F("1"),
                       self._F("x^7 + x^5 + x^4 + x^2 + 1"),
                       self._F("x^7 + x^3 + x^2 + x + 1")]
        self._sb_affine_E = sum([sb_E_coeffs[i] * self._polyring.gen()**(2**i)
                                 for i in range(8)])
        self._sb_affine_E += self._polyring("x^6 + x^5 + x + 1")
        sb_D_coeffs = [self._F("x^2 + 1"),
                       self._F("x^7 + x^6 + x^5 + x^4 + x^3 + x^2 + x"),
                       self._F("x^6 + x^5 + x^4 + x^3 + x^2 + x + 1"), 
                       self._F("x^6 + x^4 + x^3 + x"), 
                       self._F("x^6 + x^5 + x^4 + x^3"),
                       self._F("x^6 + x^4 + x^3 + 1"), 
                       self._F("x^7 + x^6 + x^4 + x^3 + x + 1"),
                       self._F("x^6 + x^5 + x^3 + x^2 + x")]
        self._sb_affine_D = sum([sb_D_coeffs[i] * self._polyring.gen()**(2**i)
                                 for i in range(8)])
        self._sb_affine_D += self._polyring("x^2 + 1")
        self._mix_col_cx = self._polyring("(x+1)*(p^3) + p^2 + p^1 + x")
        self._mix_col_dx = self._polyring("""(x^3 + x + 1)*(p^3) + 
                                          (x^3 + x^2 + 1)*(p^2) + (x^3 + 1)*p +
                                          (x^3 + x^2 + x)""")
        self._mix_col_mod = self._polyring("p^4 + 1")

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
        
    def hex_to_GF(self, H):
        r"""
        Returns the list of elements of `\GF{2^8}` corresponding to the hex
        string ``H``. Each element in `\GF{2^8}` corresponds to a unique 
        8-bit binary string, which is represented as 2 hex characters here. 
        In particular, the element 
        `a_7x^7 + a_6x^6 + a_5x^5 + a_4x^4 + a_3x^3 + a_2x^2 + a_1x^1 + a_0` 
        corresponds to the binary string '`a_7a_6a_5a_4a_3a_2a_1a_0`\`.

        INPUT:
        
        - ``H`` -- A hex string where every two characters correspond to a
          single hex value in the range `[0, 255]`.

        OUTPUT:

        - A list of elements of `\GF{2^8}` where each element corresponds to 
          the appropriate hex value in ``H``. 

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4,4)
            sage: rgf.hex_to_GF('2f')
            [x^5 + x^3 + x^2 + x + 1]
            sage: rgf.hex_to_GF('1a2b0f')
            [x^4 + x^3 + x, x^5 + x^3 + x + 1, x^3 + x^2 + x + 1]

        We can use this to build a state matrix over `\GF{2^8}` which can be
        used as input to any of the round methods: ::

            sage: state = rgf.hex_to_GF('1147659047cf663b9b0ece8dfc0bf1f0')
            sage: state = column_matrix(rgf.block_length(), 4, state)
            sage: output = rgf.shift_rows(state)
            sage: rgf.GF_to_hex(output)
            '11cfcef0470ef1909b0b653bfc47668d'

        We can use the output from this method directly if we set the 
        ``vector`` keyword to ``True`` for the round operations. ::

            sage: state = rgf.hex_to_GF('1147659047cf663b9b0ece8dfc0bf1f0') 
            sage: output = rgf.shift_rows(state, vector=True)
            sage: rgf.GF_to_hex(output) 
            '11cfcef0470ef1909b0b653bfc47668d'
        """
        if not isinstance(H, basestring) or \
           any([c not in '0123456789abcdefABCDEF' for c in H]):
            raise TypeError("keyword \'H\' must be a hex string")

        def h_to_gf(h):
            return self._F(map(int, bin(int(h, 16))[2:].zfill(8))[::-1])
        hexes = [H[2*i] + H[2*i+1] for i in range(len(H)/2)]
        return [h_to_gf(h) for h in hexes]

    def GF_to_hex(self, GF):
        r"""
        Returns the hex string representation of ``GF``.  Each element in 
        `\GF{2^8}` corresponds to a unique 8-bit binary string, which is
        represented as 2 hex characters here. In particular, the element
        `a_7x^7 + a_6x^6 + a_5x^5 + a_4x^4 + a_3x^3 + a_2x^2 + a_1x^1 + a_0`
        corresponds to the binary string '`a_7a_6a_5a_4a_3a_2a_1a_0`\`.

        INPUT:
        
        - ``GF`` -- Either a state matrix over `\GF{2^8}`, a list of elements
          from `\GF{2^8}`, or a single element from `\GF{2^8}`

        OUTPUT:

        - A hex string representation of ``GF``, where every two characters in
          the string correspond to a single element in ``\GF{2^8}``.

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
            sage: output = rgf.add_round_key(plain, key, vector=True)
            sage: rgf.GF_to_hex(output)
            'e142cd5fcd9d6d94a3340793034391b5'
        """
        from sage.rings.finite_rings.element_base import is_FiniteFieldElement
        # Test that the base field of GF is isomorphic to rgf._F
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

    def bin_to_GF(self, B):
        r"""
        Returns the list of elements of `\GF{2^8}` corresponding to the binary
        string ``B``. Each element in `\GF{2^8}` corresponds to a unique 
        8-bit binary string. In particular, the element 
        `a_7x^7 + a_6x^6 + a_5x^5 + a_4x^4 + a_3x^3 + a_2x^2 + a_1x^1 + a_0` 
        corresponds to the binary string '`a_7a_6a_5a_4a_3a_2a_1a_0`\`.

        INPUT:

        - ``B`` -- A binary string where every 8 bits correspond to a value in
          the range `[0, 255]`

        OUTPUT:

        - A list of elements of `\GF{2^8}` where each element corresponds to
          the appropriate 8-bit binary string in ``B``. 

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4,4)
            sage: bs = '01010011'
            sage: rgf.bin_to_GF()
            [x^6 + x^4 + x + 1]
            sage: bs = '000100000111001000110101110001101101011100110101'
            sage: rgf.bin_to_GF(bs)
            [x^4,
            x^6 + x^5 + x^4 + x,
            x^5 + x^4 + x^2 + 1,
            x^7 + x^6 + x^2 + x,
            x^7 + x^6 + x^4 + x^2 + x + 1,
            x^5 + x^4 + x^2 + 1]

        We can use the output of this method directly if we set the keyword
        ``vector`` to ``True`` in the round operations. ::
        
            sage: bs = '11101011100111110000000111001100' * 4
            sage: len(bs)
            128
            sage: output = rgf.sub_bytes(rgf.bin_to_GF(bs), vector=True)
            sage: rgf.GF_to_bin(output)
            '11101001110110110111110001001011111010011101101101111100010010111110100111011011011111000100101111101001110110110111110001001011'
        """
        if not isinstance(B, basestring) or \
           any([c not in '01' for c in B]):
            raise TypeError("keyword \'B\' must be a binary string")
        
        def b_to_gf(b):
            return self._F(map(int, b)[::-1])
        bins = [B[8*i : 8*(i+1)] for i in range(len(B) / 8)]
        return [b_to_gf(b) for b in bins]

    def GF_to_bin(self, GF):
        r"""
        Returns the binary string representation of ``GF``.  Each element in 
        `\GF{2^8}` corresponds to a unique 8-bit binary string. 
        In particular, the element
        `a_7x^7 + a_6x^6 + a_5x^5 + a_4x^4 + a_3x^3 + a_2x^2 + a_1x^1 + a_0`
        corresponds to the binary string '`a_7a_6a_5a_4a_3a_2a_1a_0`\`.

        INPUT:
        
        - ``GF`` -- Either a state matrix over `\GF{2^8}`, a list of elements
          from `\GF{2^8}`, or a single element from `\GF{2^8}`

        OUTPUT:

        - A binary string representation of ``GF``, where every eight 
          characters in the string correspond to a single element in
          ``\GF{2^8}``.

        EXAMPLES::

            sage: from sage.crypto.block_cipher.rijndael_gf import RijndaelGF
            sage: rgf = RijndaelGF(4,4)
            sage: F.<a> = GF(2^8)                                             
            sage: els = [a^7 + a^5 + a + 1, a^7 + a^6 + a^4 + a^2]
            sage: rgf.GF_to_bin(els) 
            '1010001111010100'

        We can use this to get clearer output from any of the round functions.
        ::

        sage: plain = '11101011100111110000000111001100' * 4
        sage: plain_state = rgf.bin_to_GF(plain)
        sage: key = '00110011100000001111100111010111' * 4
        sage: key_state = rgf.bin_to_GF(key)
        sage: output = rgf.add_round_key(plain_state, key_state, vector=True)
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

            state = column_matrix(self._Nb, 4, self.hex_to_GF(plain))
            key_state = column_matrix(self._Nk, 4, self.hex_to_GF(key))
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

            state = column_matrix(self._Nb, 4, self.bin_to_GF(plain))
            key_state = column_matrix(self._Nk, 4, self.bin_to_GF(key))
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
            if len(plain) != 8 * self._Nb:
                msg = '\'plain\' keyword\'s length must be {0}, not{1}'
                raise ValueError(msg.format(8 * self._Nb, len(ciphertext)))
            if not isinstance(key, basestring) or \
               any([c not in '0123456789abcdefABCDEF' for c in key]):
                raise TypeError("\'key\' keyword must be a hex string")
            if len(key) != 8 * self._Nk:
                msg = '\'key\' keyword\'s length must be {0}, not {1}'
                raise ValueError(msg.format(8 * self._Nk, len(key)))

            state = column_matrix(self._Nb, 4, self.hex_to_GF(ciphertext))
            key_state = column_matrix(self._Nk, 4, self.hex_to_GF(key))
            roundKeys = self.expand_key(key_state)
        elif format == 'binary':
            if not isinstance(ciphertext, basestring) or \
               any([c not in '01' for c in ciphertext]):
                raise TypeError(("\'ciphertext\' keyword must be a binary "
                                 "string"))
            if len(ciphertext) != 32 * self._Nb:
                msg = '\'plain\' keyword\'s length must be {0}, not {1}'
                raise ValueError(msg.format(32 * self._Nb, len(ciphertext)))
            if not isinstance(key, basestring) or \
               any([c not in '01' for c in key]):
                raise TypeError("\'key\' keyword must be a binary string")
            if len(key) != 32 * self._Nk:
                msg = '\'key\' keyword\'s length must be {0}, not {1}'
                raise ValueError(msg.format(32 * self._Nk, len(key)))

            state = column_matrix(self._Nb, 4, self.bin_to_GF(ciphertext))
            key_state = column_matrix(self._Nk, 4, self.bin_to_GF(key))
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

    def _check_valid_GFvector(self, GFv, keyword, expected_length=None):
        if expected_length == None:
            expected_length = self._Nb * 4
        msg = "keyword \'{0}\' must be a length {1} vector over GF({2})"
        msg = msg.format(keyword, expected_length, self._F.order())
        if not all([el.parent().order() == self._F.order() and 
                    el.parent().is_field() for el in GFv]) or \
           not len(GFv) == expected_length:
            raise TypeError(msg)

    def expand_key(self, key, vector=False):
        if vector:
            self._check_valid_GFvector(key, 'key', self._Nk*4)
        else:
            self._check_valid_GFmatrix(key, 'key', self._Nk)
        
        # Is this bad form?
        def add_cols(col1, col2):
            return map(lambda (x,y): x + y, zip(col1, col2))
        
        key_cols = []
        for i in range(self._Nb * (self._Nr + 1)):
            key_cols.append([])

        if vector:
            for j in range(self._Nk):
                for i in range(4):
                    key_cols[j].append(key[4*j + i])
        else:
            for j in range(self._Nk):
                key_cols[j] = list(key.columns()[j])

        for j in range(self._Nk, self._Nb * (self._Nr + 1)):
            if j % self._Nk == 0:
                # Apply non-linear function to k[j - 1]
                add_key = map(self._srd, key_cols[j - 1])
                add_key = self._rotate_row(add_key, 1)
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
                
    def add_round_key(self, state, round_key, vector=False):
        if vector:
            self._check_valid_GFvector(state, 'state')
            self._check_valid_GFvector(round_key, 'round_key', self._Nk*4)
            return map(lambda (x,y) : x + y, zip(state, round_key))
        else:
            self._check_valid_GFmatrix(state, 'state')
            self._check_valid_GFmatrix(round_key, 'round_key', self._Nk)
            return state + round_key

    def _srd(self, el, algorithm='encrypt'):
        if algorithm == 'encrypt':
            if el != self._F.zero():
                newEl = el ** -1
            else:
                newEl = el

            return self._sb_affine_E.subs(p = newEl)
        elif algorithm == 'decrypt':
            newEl = self._sb_affine_D.subs(p = el)

            if newEl == self._F.zero():
                return newEl
            else:
                return newEl ** -1
        else:
            raise ValueError(("keyword 'algorithm' must be either 'encrypt' "
                             "or 'decrypt'"))
    
    def sub_bytes(self, state, algorithm='encrypt', vector=False):
        if vector:
            self._check_valid_GFvector(state, 'state')
            return [self._srd(el, algorithm) for el in state]
        else:
            self._check_valid_GFmatrix(state, 'state')

            new_columns = []
            for col in state.columns():
                new_columns.append([self._srd(el, algorithm) for el in col])
            return column_matrix(new_columns)

    def mix_columns(self, state, algorithm='encrypt', vector=False):
        if vector:
            self._check_valid_GFvector(state, 'state')
        else:
            self._check_valid_GFmatrix(state, 'state')

        if algorithm == 'encrypt':
            constant = self._mix_col_cx
        elif algorithm == 'decrypt':
            constant = self._mix_col_dx
        else:
            raise ValueError(("keyword 'algorithm' must be either 'encrypt' "
                             "or 'decrypt'"))
        if vector:
            newState = []
            for i in range(len(state)/4):
                bx = [state[j*4:j*4 + i] for i in range(4)]
                bx = constant * self._polyring(list(col))
                bx = bx + [self._F.zero()]*(4 - len(bx))
                newState += bx
            return newState
        else:
            newState = [[], [], [], []]
            new_columns = []
            for col in state.columns():
                bx = constant * self._polyring(list(col))
                bx = list(bx.mod(self._mix_col_mod))
                bx = bx + [self._F.zero()]*(4 - len(bx))
                new_columns.append(bx)
            return column_matrix(new_columns)

    def _rotate_row(self, row, n):
        row = list(row)
        return row[n:] + row[:n]

    def shift_rows(self, state, algorithm='encrypt', vector=False):
        if vector:
            self._check_valid_GFvector(state, 'state')
            state = column_matrix(self._Nb, 4, state)   
        else:
            self._check_valid_GFmatrix(state, 'state')

        if algorithm == 'encrypt':
            offsets = self._shiftrows_offsets_E
        elif algorithm == 'decrypt':
            offsets = self._shiftrows_offsets_D
        else:
            raise ValueError(("keyword 'algorithm' must be either 'encrypt' "
                             "or 'decrypt'"))

        rows = [[], [], [], []]
        rows[0] = self._rotate_row(state[0], offsets[self._Nb - 4][0])
        rows[1] = self._rotate_row(state[1], offsets[self._Nb - 4][1])
        rows[2] = self._rotate_row(state[2], offsets[self._Nb - 4][2])
        rows[3] = self._rotate_row(state[3], offsets[self._Nb - 4][3])
        if vector:
            vec = []
            for col in matrix(rows).columns():
                vec += list(col)
            return vec
        else:
            return matrix(rows)
