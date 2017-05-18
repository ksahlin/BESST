'''
    Created on Jun 2, 2012
    
    @author: ksahlin
    This module calculates the expected number of links spanning two contigs 
    of length len1 and len2 with a gap of length d. This formula was 
    first published by Zerbino et al. "Pebble and Rock Band: Heuristic 
    Resolution of Repeats and Scaffolding in the Velvet Short-Read de Novo Assembler",
    PlosOne 2009. 
    This module adds parameters for read length and softclipped bases,
    See info in Param-class.
'''

import argparse
import math
from decimal import Decimal, getcontext

def normcdf(x, mu, sigma):
    t = x - mu;
    y = 0.5 * math.erfc(-t / (sigma * math.sqrt(2.0)));
    if y > 1.0:
        y = 1.0;
    return y

def normpdf(x, mu, sigma):
    #Get much better approximations with Decimal (simply more decimals)
    getcontext().prec = 100
    u = Decimal(str(x - mu)) / Decimal(str(abs(sigma)))
    y = float(str((1 / Decimal(str((math.sqrt(2 * math.pi) * abs(sigma))))) * Decimal(str(-u * u / 2)).exp()))
    return y



class Param(object):
    """Holds parameter values neccessary to calculate expected number of
    links between two contigs. Parameters that are the same over a whole
    assembly run are contained here. Parameters that vary for each problem 
    instance are sent explicitly to the ExpectedLinks() function. 
    Attributes:
    mean            -- mean fragment length of the library
    stddev          -- Standard deviation of the library
    read_len        -- The mate length (i.e. the length of the sequences part)
                        Assuming same length of both mates
    cov             -- The mean coverage for the library
    softclipped     -- Number of allowed bases from read "hanging" outside contig.
    readfrequency   -- Expected number of positions between two consecutive PE/MP-fragments.
                        See example below
    
    Readfrequency is the expected number of positions from which a new 
    read starts. i.e. readfrequency = 2*read_len/mean_coverage   
    (read_len is the nr of nucleotides sampled from each side of a fragment,
    assume equally long on each side.
    example: read_len = 100, coverage=100 -> k = 2, 
    i.e a new read starts on average every 2nd position on the genome.

    (read_len - softclipped) ultimately means the number ofbases from a read that needs to reside
    within a contig in order to be considered as mapped to that contig.
    """
    def __init__(self, mean, stddev, cov, read_len, softclipped):
        self.mean = mean
        self.stddev = stddev
        self.read_len = read_len
        self.cov = cov
        self.softclipped = softclipped
        self.readfrequency = 2 * self.read_len / self.cov


def ExpectedLinks(len1, len2, d, param):
    '''
    Caculates the expected number of links spanning two contigs of length 
    len1 and len2 with a gap of length d. 
    '''
    std_dev = float(param.stddev)
    # We should not expect more links for a negative gap since the negative gap can come from misassembly,
    # bad gap estimate, overlapping k-mers from cutting assembly graph or haplotypic regions 
    # Even if the overlap is correct we should not expect too many links in the overlapping region since part of
    # the PE/MP will map to the same contig so this will create a problem anyway, 
    # therefore we count with gap = max(d,0)
    gap = max(d, 0)
    #Specifying input arguments
    b1 = (len1 + len2 + gap + 2 * param.softclipped - param.mean) / std_dev
    a1 = (max(len1, len2) + gap + (param.read_len - param.softclipped) - param.mean) / std_dev
    b2 = (min(len1, len2) + gap + (param.read_len - param.softclipped) - param.mean) / std_dev
    a2 = (gap + 2 * (param.read_len - param.softclipped) - param.mean) / std_dev
    def Part(a, b):
        expr1 = (min(len1, len2) - (param.read_len - param.softclipped)) / param.readfrequency * normcdf(a, 0, 1)
        expr2 = -(-param.softclipped) / param.readfrequency * normcdf(b, 0, 1)
        expr3 = (b * std_dev) / param.readfrequency * (normcdf(b, 0, 1) - normcdf(a, 0, 1))
        expr4 = (std_dev / param.readfrequency) * (normpdf(b, 0, 1) - normpdf(a, 0, 1))
        value = expr1 + expr2 + expr3 + expr4
        return value
    E_links = Part(a1, b1) - Part(a2, b2)
    return(E_links)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser(description="calculates the expected number of links spanning two contigs. ")
    arg_parser.add_argument("mean", type=int, help="Mean insert size.")
    arg_parser.add_argument("stddev", type=int, help="Standard deviation of insert size.")
    arg_parser.add_argument("cov", type=int, help="Mean coverage.")
    arg_parser.add_argument("readlen", type=int, help="Read length.")
    arg_parser.add_argument("soft", type=int, help="Number of softclipped bases allowed.")
    arg_parser.add_argument("len1", type=int, help="Contig1 length.")
    arg_parser.add_argument("len2", type=int, help="Contig2 length.")
    arg_parser.add_argument("d", type=int, help="Gap size")



    args = arg_parser.parse_args()

    param = Param(args.mean, args.stddev, args.cov, args.readlen, args.soft)
    print ExpectedLinks(args.len1, args.len2, args.d, param)