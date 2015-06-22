r"""
Rijndael-GF

A strictly algebraic implementation of the Rijndael Cipher (more commonly 
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
        self._shiftrows_offset_E = matrix([[0,1,2,3], [0,1,2,3], [0,1,2,3],
                                   [0,1,2,4], [0,1,3,4]])
        self._shiftrows_offset_D = matrix([[0,-1,-2,-3], [0,-1,-2,-3],
                                           [0,-1,-2,-3], [0,-1,-2,-4],
                                           [0,-1,-3,-4]])
        self._Nr = self._round_num_table[self._Nb - 4][self._Nk - 4]
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
        msg = """Rijndael-GF block cipher with block length {0}, key length 
              {1}, and {2} rounds."""
        return msg.format(self._Nb, self._Nk, self._Nr)

    def block_length(self):
        return self._Nb

    def key_length(self):
        return self._Nk

    def number_rounds(self):
        return self._Nk
        
    def _hex_to_GF(self, h):
        return self._F(map(int, bin(int(h,16))[2:].zfill(8))[::-1])

    def _GF_to_hex(self, gf):
        return hex(gf.integer_representation())[2:].zfill(2)

    def _bin_to_GF(self, b):
        return self._F(map(int, b)[::-1])

    def _GF_to_bin(self, gf):
        return bin(gf.integer_representation())[2:].zfill(8)

    def _binString_to_hexString(self, binS):
        btw = self._bits_to_word
        bin_vals = []
        for i in range(len(binS) / btw):
            bin_vals.append(''.join([binS[i*btw + j] for j in range(btw)]))
        return ''.join([hex(int(b, 2))[2:].zfill(2) for b in bin_vals])

    def _hex_to_Fmatrix(self, hexS):
        hexes = [hexS[2*i] + hexS[2*i+1] for i in range(len(hexS)/2)]
        columns = []

        for i in range(self._Nb):
            columns.append([self._hex_to_GF(hexes[i*4 + j]) for j in range(4)])
        return column_matrix(columns)

    def _Fmatrix_to_hex(self, fmatrix):
        return ''.join([self._GF_to_hex(el) 
                        for col in fmatrix.columns() for el in col])

    def _bin_to_Fmatrix(self, binS):
        return self._hex_to_Fmatrix(self._binString_to_hexString(binS))

    def _Fmatrix_to_bin(self, fmatrix):
        return ''.join([self._GF_to_bin(el)
                        for col in fmatrix.columns() for el in col])
        
    def encrypt(self, plain, key, format='hex'):
        if format == 'hex':
            if not isinstance(plain, basestring) or \
               any([c not in '0123456789abcdefABCDEF' for c in plain]):
                raise ValueError("\'plain\' keyword must be a hex string")
            if len(plain) != 8 * self._Nb:
                msg = '\'plain\' keyword\'s length must be {0}, not{1}'
                raise ValueError(msg.format(8 * self._Nb, len(plain)))
            if not isinstance(key, basestring) or \
               any([c not in '0123456789abcdefABCDEF' for c in key]):
                raise ValueError("\'key\' keyword must be a hex string")
            if len(key) != 8 * self._Nk:
                msg = '\'key\' keyword\'s length must be {0}, not {1}'
                raise ValueError(msg.format(8 * self._Nk, len(key)))
            state = self._hex_to_Fmatrix(plain)
            roundKeys = self.expand_key(self._hex_to_Fmatrix(key))
        elif format == 'binary':
            if len(plain) != 32 * self._Nb:
                msg = '\'plain\' keyword\'s length must be {0}, not {1}'
                raise ValueError(msg.format(32 * self._Nb, len(plain)))
            if not isinstance(plain, basestring) or \
               any([c not in '01' for c in plain]):
                raise ValueError("\'plain\' keyword must be a binary string")
            if len(key) != 32 * self._Nk:
                msg = '\'key\' keyword\'s length must be {0}, not {1}'
                raise ValueError(msg.format(32 * self._Nk, len(key)))
            if not isinstance(key, basestring) or \
               any([c not in '01' for c in key]):
                raise ValueError("\'key\' keyword must be a binary string")
            state = self._bin_to_Fmatrix(plain)
            roundKeys = self.expand_key(self._bin_to_Fmatrix(key))
        else:
            raise ValueError("""\'format\' keyword must be either \'hex\' or
                             \'binary\'""")

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
            return self._Fmatrix_to_hex(state)
        else:
            return self._Fmatrix_to_bin(state)

    def decrypt(self, ciphertext, key, format='hex'):
        if format == 'hex':
            if not isinstance(ciphertext, basestring) or \
               any([c not in '0123456789abcdefABCDEF' for c in ciphertext]):
                raise ValueError("\'ciphertext\' keyword must be a hex string")
            if len(plain) != 8 * self._Nb:
                msg = '\'plain\' keyword\'s length must be {0}, not{1}'
                raise ValueError(msg.format(8 * self._Nb, len(ciphertext)))
            if not isinstance(key, basestring) or \
               any([c not in '0123456789abcdefABCDEF' for c in key]):
                raise ValueError("\'key\' keyword must be a hex string")
            if len(key) != 8 * self._Nk:
                msg = '\'key\' keyword\'s length must be {0}, not {1}'
                raise ValueError(msg.format(8 * self._Nk, len(key)))
            state = self._hex_to_Fmatrix(plain)
            roundKeys = self.expand_key(self._hex_to_Fmatrix(key))
        elif format == 'binary':
            if len(ciphertext) != 32 * self._Nb:
                msg = '\'plain\' keyword\'s length must be {0}, not {1}'
                raise ValueError(msg.format(32 * self._Nb, len(ciphertext)))
            if not isinstance(ciphertext, basestring) or \
               any([c not in '01' for c in ciphertext]):
                raise ValueError("""\'ciphertext\' keyword must be a binary
                                 string""")
            if len(key) != 32 * self._Nk:
                msg = '\'key\' keyword\'s length must be {0}, not {1}'
                raise ValueError(msg.format(32 * self._Nk, len(key)))
            if not isinstance(key, basestring) or \
               any([c not in '01' for c in key]):
                raise ValueError("\'key\' keyword must be a binary string")
            state = self._bin_to_Fmatrix(ciphertext)
            roundKeys = self.expand_key(self._bin_to_Fmatrix(key))
        else:
            raise ValueError("""\'format\' keyword must be either \'hex\' or
                             \'binary\'""")

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
            return self._Fmatrix_to_hex(state)
        else:
            return self._Fmatrix_to_bin(state)
                  
    def _check_valid_GFmatrix(self, GFm, keyword, expected_cols = self._Nb):
        msg = "keyword \'{0}\' must be a {1} x {2} matrix over GF({3})"
        msg = msg.format(keyword, 4, expected_cols, self._F.order())
        if not isinstance(GFm, Matrix_dense) or \
           not (GFm.base_ring().order() == self._F.order() and \
                GFm.base_ring().is_field()) or \
           not (GFm.ncols() == expected_cols and GFm.nrows() == 4):
            raise TypeError(msg)

    def expand_key(self, keyMatrix):
        self._check_valid_GFmatrix(keyMatrix, 'keyMatrix', self._Nk)
        
        # Is this bad form?
        def add_cols(col1, col2):
            return map(lambda (x,y): x + y, zip(col1, col2))
        
        key_cols = []
        for i in range(self._Nb * (self._Nr + 1)):
            key_cols.append([])
        
        for j in range(self._Nk):
            key_cols[j] = list(keyMatrix.columns()[j])

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
                

    def add_round_key(self, state, round_key):
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
            raise ValueError("""keyword 'algorithm' must be either 'encrypt'
                             or 'decrypt'""")
    
    def sub_bytes(self, state, algorithm='encrypt'):
        self._check_valid_GFmatrix(state, 'state')

        new_columns = []
        for col in state.columns():
            new_columns.append([self._srd(el, algorithm) for el in col])
        return column_matrix(new_columns)

    def mix_columns(self, state, algorithm='encrypt'):
        self._check_valid_GFmatrix(state, 'state')
        if algorithm == 'encrypt':
            constant = self._mix_col_cx
        elif algorithm == 'decrypt':
            constant = self._mix_col_dx
        else:
            raise ValueError("""keyword 'algorithm' must be either 'encrypt'
                             or 'decrypt'""")

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

    def shift_rows(self, state, algorithm='encrypt'):
        self._check_valid_GFmatrix(state, 'state')
        if algorithm == 'encrypt':
            offsets = self._shiftrows_offsets_E
        elif algorithm == 'decrypt':
            offsets = self._shiftrows_offsets_D
        else:
            raise ValueError("""keyword 'algorithm' must be either 'encrypt'
                             or 'decrypt'""")

        rows = [[], [], [], []]
        rows[0] = self._rotate_row(state[0], offsets[self._Nb - 4][0])
        rows[1] = self._rotate_row(state[1], offsets[self._Nb - 4][1])
        rows[2] = self._rotate_row(state[2], offsets[self._Nb - 4][2])
        rows[3] = self._rotate_row(state[3], offsets[self._Nb - 4][3])
        return matrix(rows)
