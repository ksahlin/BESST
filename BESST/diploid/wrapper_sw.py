'''
Created on Jun 26, 2013

@author: ksahlin
'''

import ctypes
class ReturnValues(ctypes.Structure):
    _fields_ = [("score", ctypes.c_int),
                ("i_max", ctypes.c_int),
                ("j_max", ctypes.c_int)]

def wrap_sw(seq1, seq2, mu, delta):
    c_s1 = ctypes.c_char_p(seq1)
    c_s2 = ctypes.c_char_p(seq2)
    res = ReturnValues(0, 0, 0)
    SW_fcn = ctypes.CDLL('BESST/diploid/lib/smith_waterman.so') #BESST/diploid/
    SW_fcn.SW(c_s1, c_s2, ctypes.byref(res)) # linux or when mingw used on windows
    return(res.score, res.i_max, res.j_max)

if __name__ == '__main__':
    wrap_sw('AGAAGGAAGGGGGAGAGTTTG',
            'AGAAGGAAGGTTGGGAGAGTTTG' , 0, 0)
