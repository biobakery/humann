"""
HUMAnN2: chi2cdf module
Computes the chi-square cumulative distribution

Copyright (c) 2014 Harvard School of Public Health

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""

import math
import re
import sys

# Adapted from samtools; will be occasionally inaccurate due to iteration stoppage
def incomplete_gamma1( dS, dZ ):
    
    dS, dZ = (float(d) for d in (dS, dZ))
    dSum = dX = 1
    for i in range( 1, 10000 ):
        dX *= dZ / ( dS + i )
        dSum += dX
        if ( dX / dSum ) < 1e-14:
            break
    return ( math.exp( ( dS * math.log( dZ ) ) - dZ - _log_gamma( dS + 1 ) + math.log( dSum ) ) if dZ else 0 )

def _log_gamma( dZ ):
    
    dX = 0
    dX += 0.1659470187408462e-06 / ( dZ + 7 )
    dX += 0.9934937113930748e-05 / ( dZ + 6 )
    dX -= 0.1385710331296526     / ( dZ + 5 )
    dX += 12.50734324009056      / ( dZ + 4 )
    dX -= 176.6150291498386      / ( dZ + 3 )
    dX += 771.3234287757674      / ( dZ + 2 )
    dX -= 1259.139216722289      / ( dZ + 1 )
    dX += 676.5203681218835      / dZ
    dX += 0.9999999999995183;
    return ( math.log( dX ) - 5.58106146679532777 - dZ + ( ( dZ - 0.5 ) * math.log( dZ + 6.5 ) ) )

# Implementation thanks to http://www.crbond.com/math.htm thanks to Zhang and Jin
# Modified to only return normalized/regularized lower incomplete gamma
def incomplete_gamma2( dA, dX ):

    if ( dA < 0 ) or ( dX < 0 ):
        return None
    if not dX:
        return 0
    xam = -dX + dA * math.log( dX )
    if ( xam > 700 ) or ( dA > 170 ):
        return 1
    if dX <= ( dA + 1 ):
        r = s = 1.0 / dA
        for k in range( 1, 61 ):
            r *= float(dX) / ( dA + k )
            s += r
            if abs( r / s ) < 1e-15:
                break
        ga = math.gamma( dA )
        gin = math.exp( xam ) * s
        return ( gin / ga )

    t0 = 0
    for k in range( 60, 0, -1 ):
        t0 = float(k - dA) / ( 1 + ( float(k) / ( dX + t0 ) ) )
    gim = math.exp( xam ) / ( dX + t0 )
    ga = math.gamma( dA )
    return ( 1 - ( gim / ga ) )

def chi2cdf( dX, dK ):
    
    dK, dX = (( d / 2 ) for d in (dK, dX))
    dRet = incomplete_gamma1( dK, dX )
    if abs( dRet ) != float("Inf"):
        return dRet
    return incomplete_gamma2( dK, dX )

