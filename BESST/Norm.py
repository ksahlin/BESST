'''
    Created on Jun 2, 2012
    
    @author: ksahlin
    
    This file is part of BESST.
    
    BESST is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    BESST is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with BESST.  If not, see <http://www.gnu.org/licenses/>.
    '''


from math import exp, sqrt, pi
from decimal import *
def erf(x):
    ## error function approximation with an error less than 1.5 * 10-7 for all x
    # save the sign of x
    sign = 1 if x >= 0 else -1
    x = abs(x)

    # constants
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    p = 0.3275911

    # A&S formula 7.1.26
    t = 1.0 / (1.0 + p * x)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * exp(-x * x)
    return sign * y # erf(-x) = -erf(x)

def erfcc(x):
    """Complementary error function."""
    z = abs(x)
    t = 1. / (1. + 0.5 * z)
    r = t * exp(-z * z - 1.26551223 + t * (1.00002368 + t * (.37409196 +
        t * (.09678418 + t * (-.18628806 + t * (.27886807 +
        t * (-1.13520398 + t * (1.48851587 + t * (-.82215223 +
        t * .17087277)))))))))
    if (x >= 0.):
        return r
    else:
        return 2. - r

def normcdf(x, mu, sigma):
    t = x - mu;
    y = 0.5 * erfcc(-t / (sigma * sqrt(2.0)));
    if y > 1.0:
        y = 1.0;
    return y

def normpdf(x, mu, sigma):
    #Get much better approximations with Decimal (simply more decimals)
    getcontext().prec = 100
    u = Decimal(str(x - mu)) / Decimal(str(abs(sigma)))
    y = float(str((1 / Decimal(str((sqrt(2 * pi) * abs(sigma))))) * Decimal(str(-u * u / 2)).exp()))
    return y
