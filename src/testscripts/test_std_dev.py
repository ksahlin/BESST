import sys
from math import exp,sqrt,pi
from decimal import *

def erf(x):
    ## error function approximation with an error less than 1.5 * 10-7 for all x
    # save the sign of x
    sign = 1 if x >= 0 else -1
    x = abs(x)

    # constants
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911

    # A&S formula 7.1.26
    t = 1.0/(1.0 + p*x)
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x)
    return sign*y # erf(-x) = -erf(x)

def erfcc(x):
    """Complementary error function."""
    z = abs(x)
    t = 1. / (1. + 0.5*z)
    r = t * exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+
        t*(.09678418+t*(-.18628806+t*(.27886807+
        t*(-1.13520398+t*(1.48851587+t*(-.82215223+
        t*.17087277)))))))))
    if (x >= 0.):
        return r
    else:
        return 2. - r
    
#def normcdf(x, mu, sigma):
#    t = x-mu;
#    y = 0.5*erfcc(-t/(sigma*sqrt(2.0)));
#    if y>1.0:
#        y = 1.0;
#    return y
#
def normpdf(x, mu, sigma):
    getcontext().prec = 100
    u = Decimal((x-mu))/Decimal(abs(sigma))
    y = float(str((1/Decimal(str((sqrt(2*pi)*abs(sigma)))))*Decimal(str(-u*u/2)).exp()))
    #y = (1/(sqrt(2*pi)*abs(sigma)))*exp(-u*u/2)
    return y



def E_O(d,mean,stdDev,c1Len,c2Len,readLen):

    def CalcG_prime_d(d,c_min,c_max,c1Len,c2Len,readLen):
        
        term1 = - 0.5 * ( erf((c_min+c_max+d+1-mean)/(2**0.5*float(stdDev))) + erf((d +2*readLen - 1 - mean)/(2**0.5*float(stdDev))) )
        term2 = + 0.5 * ( erf((d + c_max + readLen - mean)/(2**0.5*float(stdDev))) + erf((d + c_min + readLen - mean)/(2**0.5*float(stdDev)))  )
        g_prime_d = term1 + term2
        return g_prime_d
    
    def CalcGd(d,c1Len,c2Len,readLen):
        #term 1,2 and 3 denodes what part of the function we are integrating term1 for first (ascending), etc...
        term1 = (c_min-readLen+1)/2.0* (erf((c_max+d+readLen-mean)/((2**0.5)*stdDev))- erf((c_min+d+readLen-mean)/((2**0.5)*stdDev))   )
        
        term2 = (c_min+c_max+d-mean+1)/2.0 *( erf((c_min+c_max+d+1-mean)/(2**0.5*float(stdDev))) - erf((c_max+readLen+d-mean)/(2**0.5*float(stdDev))) )
        
        term3 = (d+2*readLen-mean-1)/2.0 *( erf((d+2*readLen-1-mean)/(2**0.5*float(stdDev))) - erf((c_min+d+readLen-mean)/((2**0.5)*stdDev)) )
        
        term4 = stdDev / ((2*pi)**0.5) * ( 2.718**(-( (c_min+c_max+d+1-mean)**2)/(float(2*stdDev**2))) + 2.718**(-( (d+2*readLen-1-mean)**2)/(float(2*stdDev**2))) )
        
        term5 = - stdDev / ((2*pi)**0.5) * ( 2.718**(-( (c_max+readLen+d-mean)**2)/(float(2*stdDev**2))) + 2.718**(-( (c_min+readLen+d-mean)**2)/(float(2*stdDev**2))) )
        
        g_d = term1 + term2 + term3+ term4 + term5
        return g_d
    
    c_min=min(c1Len,c2Len)
    c_max=max(c1Len,c2Len)
    g_prime_d = CalcG_prime_d(d,c_min,c_max,c1Len,c2Len,readLen)
    g_d = CalcGd(d,c1Len,c2Len,readLen)
    a = stdDev**2*g_prime_d + mean*g_d   
    E_o = a / g_d  - d
    #print 'E[O|d ='+str(d)+'] =', E_o
    #print 'Theoretical base case mean:',  mean + stdDev**2/ float(mean + 1) 
    #print 'Theoretical base case std_dev:', ( stdDev**2 - stdDev**4/ (mean + 1)**2 )**0.5
    return E_o

def E_O_square(d,mean,stdDev,c1Len,c2Len,readLen):

    def CalcG_prime_d(d,c_min,c_max,c1Len,c2Len,readLen):
        
        term1 = - 0.5 * ( erf((c_min+c_max+d+1-mean)/(2**0.5*float(stdDev))) + erf((d +2*readLen - 1 - mean)/(2**0.5*float(stdDev))) )
        term2 = + 0.5 * ( erf((d + c_max + readLen - mean)/(2**0.5*float(stdDev))) + erf((d + c_min + readLen - mean)/(2**0.5*float(stdDev)))  )
        g_prime_d = term1 + term2
        return g_prime_d
    
    def CalcGd(d,c1Len,c2Len,readLen):
        #term 1,2 and 3 denodes what part of the function we are integrating term1 for first (ascending), etc...
        term1 = (c_min-readLen+1)/2.0* (erf((c_max+d+readLen-mean)/((2**0.5)*stdDev))- erf((c_min+d+readLen-mean)/((2**0.5)*stdDev))   )
        
        term2 = (c_min+c_max+d-mean+1)/2.0 *( erf((c_min+c_max+d+1-mean)/(2**0.5*float(stdDev))) - erf((c_max+readLen+d-mean)/(2**0.5*float(stdDev))) )
        
        term3 = (d+2*readLen-mean-1)/2.0 *( erf((d+2*readLen-1-mean)/(2**0.5*float(stdDev))) - erf((c_min+d+readLen-mean)/((2**0.5)*stdDev)) )
        
        term4 = stdDev / ((2*pi)**0.5) * ( 2.718**(-( (c_min+c_max+d+1-mean)**2)/(float(2*stdDev**2))) + 2.718**(-( (d+2*readLen-1-mean)**2)/(float(2*stdDev**2))) )
        
        term5 = - stdDev / ((2*pi)**0.5) * ( 2.718**(-( (c_max+readLen+d-mean)**2)/(float(2*stdDev**2))) + 2.718**(-( (c_min+readLen+d-mean)**2)/(float(2*stdDev**2))) )
        g_d = term1 + term2 + term3+ term4 + term5
        return g_d
    
    def CalcB(d,c_min,c_max,c1Len,c2Len,readLen,g_d,g_prime_d):
        c1 = 0 # norm.pdf( c_min + c_max + d + 1 , mean, stdDev) + norm.pdf( d + 2*readLen - 1, mean, stdDev) 
        c2 = normpdf( d + 2*readLen - 1, mean, stdDev) - normpdf( c_min + d + readLen , mean, stdDev) - normpdf( c_max + d + readLen , mean, stdDev) + normpdf( c_min + c_max + d + 1 , mean, stdDev) 
        
        b = stdDev**4 *(c1 + c2) + mean**2*g_d + stdDev**2*g_d  + 2*stdDev**2*mean*g_prime_d 
        #print g_prime_d,g_d, c1, c2, stdDev**4 *(c1 + c2)/g_d,  b,stdDev**4 *(c1 + c2)
        return b,c1,c2
    
    c_min=min(c1Len,c2Len)
    c_max=max(c1Len,c2Len)
    g_prime_d = CalcG_prime_d(d,c_min,c_max,c1Len,c2Len,readLen)
    g_d = CalcGd(d,c1Len,c2Len,readLen)
    a = stdDev**2*g_prime_d + mean*g_d   
    b,c1,c2 = CalcB(d,c_min,c_max,c1Len,c2Len,readLen,g_d,g_prime_d)
    E_o_square = b / g_d - 2*d*(stdDev**2*(g_prime_d/g_d) + mean) + d**2 
    #print 'E[O^2|d ='+str(d)+'] =', E_o_square #, (c1+c2)*stdDev**4/g_d, stdDev**4*(g_prime_d/g_d)**2#,- 2*d*(stdDev**2*(g_prime_d/g_d) + mean), (g_prime_d/g_d + mean)  , b / g_d , d**2
    
    #print 'Shortcut var:', stdDev**2 - stdDev**4*(g_prime_d/g_d)**2 + (c1+c2)*stdDev**4/g_d
    #print stdDev**4, (g_prime_d/g_d)**2, (c1+c2)
    #print 'Theoretical base case mean:',  mean + stdDev**2/ float(mean + 1) 
    #print 'Theoretical base case std_dev:', ( stdDev**2 - stdDev**4/ (mean + 1)**2 )**0.5
    return E_o_square



if __name__ == "__main__":
    mean,stdDev,readLen,c1Len,c2Len = 1500,500,50,3000,3000
    for k in range(0,3000,500):
        e_o = E_O(k,mean,stdDev,c1Len,c2Len,readLen)
        e_o_square = E_O_square(k,mean,stdDev,c1Len,c2Len,readLen)
        #print e_o**2
        print 'Std_dev(O|d ='+str(k)+')=', (e_o_square - e_o**2)**0.5
    
    
    
    
    