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
        z = polygen(Integers(2))
        mod = z**8 + z**4 + z**3 + z + 1
        self._F = FiniteField(2**self._bits_to_word, 'x', modulus=mod)
        # Note: Any matrix created with these must be transposed!
        # TODO: See if we can have MatrixSpace read a list column-wise
        self._state_ms = MatrixSpace(self._F, 4, self._Nb)
        self._key_ms = MatrixSpace(self._F, 4, self._Nk)
        self._round_num_table = matrix([[10,11,12,13,14], [11,11,12,13,14],
                                       [12,12,12,13,14], [13,13,13,13,14],
                                       [14,14,14,14,14]])
        self._sr_offset = matrix([[0,1,2,3], [0,1,2,3], [0,1,2,3],
                                   [0,1,2,4], [0,1,3,4]])
        self._Nr = self._round_num_table[self._Nb - 4][self._Nk - 4]
        self._polyring = PolynomialRing(self._F, 'p')
        coeffs = [self._F("x^2 + 1"), 
                  self._F("x^3 + 1"), 
                  self._F("x^7 + x^6 + x^5 + x^4 + x^3 + 1"), 
                  self._F("x^5 + x^2 + 1"), 
                  self._F("x^7 + x^6 + x^5 + x^4 + x^2"), 
                  self._F("1"),
                  self._F("x^7 + x^5 + x^4 + x^2 + 1"),
                  self._F("x^7 + x^3 + x^2 + x + 1")]
        # Fix this line!
        self._sb_affine_E = sum([coeffs[i] * self._polyring.gen()**(2**i) for i in range(8)])
        self._sb_affine_E += self._polyring("x^6 + x^5 + x + 1")
        self._mix_col_cx = self._polyring("(x+1)*(p^3) + p^2 + p^1 + x")
#        self._mix_col_dx =  
        self._mix_col_mod = self._polyring("p^4 + 1")
        
    # Is there any way to make this prettier?
    def _hex_to_Fmatrix(self, hexS):
        hexes = [hexS[2*i] + hexS[2*i+1] for i in range(len(hexS)/2)]
        rows = [[], [], [], []]
        # We use modcounter here to read entries in by column, not by rows
        mod_counter = 0
        for h in hexes:
            Fval = self._F(map(int, bin(int(h,16))[2:].zfill(8))[::-1])
            rows[mod_counter].append(Fval)
            mod_counter = (mod_counter + 1) % 4

        if len(hexes)/4 == self._Nb:
            state_matrix = self._state_ms(rows)
        else:
            state_matrix = self._key_ms(rows)
        return state_matrix

    def _Fmatrix_to_hex(self, fmatrix):
        # Are inline functions like this bad form?
        def fToHex(el):
            return hex(el.integer_representation())[2:].zfill(2)
        hexs = ''.join([fToHex(el) for col in fmatrix.columns() for el in col])
        return hexs
        
    def encrypt(self, plain, key, format='hex'):
        if format == 'hex':
            if len(plain) % 2 == 1:
                raise ValueError('\'plain\' keyword''s length must be even')
            if len(key) % 2 == 1:
                raise ValueError('\'key\' keyword\'s length must be even')
            
            state = self._hex_to_Fmatrix(plain)
            roundKeys = self.expand_key(self._hex_to_Fmatrix(key))
            
        state = self.add_round_key(state, roundKeys[0])

        for r in range(self._Nr-1):
            state = self.sub_bytes(state)
            state = self.shift_rows(state)
            state = self.mix_columns(state)
            state = self.add_round_key(state, roundKeys[r+1])
            
        # Final round
        state = self.sub_bytes(state)
        state = self.shift_rows(state)
        state = self.add_round_key(state, roundKeys[self._Nr])
            
        # convert back to form given and return
        return self._Fmatrix_to_hex(state)
                  
    # We should expand a key into a list-matrix of the appropriate size,
    # then break it into 4xNb sized matrices. Return a list of these 4xNb
    # sized matrices, one for each round.
    # TODO: Make this method not so god-forsaken ugly!!!
    # TODO: Make this method work properly for Nk > 6
    def expand_key(self, keyMatrix):
        expanded = []
        for i in range(4):
            expanded.append([self._F.zero()]*(self._Nb * (self._Nr + 1)))
    
        for j in range(self._Nk):
            for i in range(4):
                expanded[i][j] = keyMatrix[i][j]
            
        for j in range(self._Nk, self._Nb * (self._Nr + 1)):
            if j % self._Nk == 0:
                expanded[0][j] = expanded[0][j - self._Nk] + self._srd(expanded[1][j-1]) + \
                                 self._F.gen()**(int(j/self._Nk)-1)
                for i in range(1,4):
                    expanded[i][j] = expanded[i][j - self._Nk] + \
                                     self._srd(expanded[(i+1) % 4][j-1])
            else:
                for i in range(4):
                    expanded[i][j] = expanded[i][j - self._Nk] + expanded[i][j-1]

        round_keys = []
        for r in range(self._Nr + 1):
            rk_entries = []
            for i in range(4):
                for j in range(self._Nb):
                    rk_entries.append(expanded[i][(r * self._Nb) + j])
            round_keys.append(self._state_ms(rk_entries))
        
        return round_keys

    def add_round_key(self, state, round_key):
        return state + round_key

    def _srd(self, el, algorithm='encrypt'):
        if algorithm == 'encrypt':
            if el != self._F.zero():
                newEl = el ** -1
            else:
                newEl = el

            return self._sb_affine_E.subs(p = newEl)
        else:
            newEl = self._sb_affine_D.subs(p = el)

            if newEl == F.zero():
                return newEl
            else:
                return newEl ** -1
    
    def sub_bytes(self, state, algorithm='encrypt'):
        new_state = []
        for row in state:
            new_state.append([self._srd(el, algorithm) for el in row])
        return self._state_ms(new_state)

    def mix_columns(self, state, algorithm='encrypt'):
        newState = [[], [], [], []]
        for col in state.columns():
            bx = self._mix_col_cx * self._polyring(list(col))
            bx = list(bx.mod(self._mix_col_mod))
            bx = bx + [self._F.zero()]*(4 - len(bx))
            for i in range(4):
                newState[i].append(bx[i])
        return self._state_ms(newState)

    def _rotate_row(self, row, n):
        row = list(row)
        return row[n:] + row[:n]

    def shift_rows(self, state, algorithm='encrypt'):
        rows = [[], [], [], []]
        rows[0] = self._rotate_row(state[0], self._sr_offset[self._Nb - 4][0])
        rows[1] = self._rotate_row(state[1], self._sr_offset[self._Nb - 4][1])
        rows[2] = self._rotate_row(state[2], self._sr_offset[self._Nb - 4][2])
        rows[3] = self._rotate_row(state[3], self._sr_offset[self._Nb - 4][3])
        return self._state_ms(rows)

# TODO: add self._srd_affine_D, 
# self.decrypt, more conversion methods, look at changing word length
