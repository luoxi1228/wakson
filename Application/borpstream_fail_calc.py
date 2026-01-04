#!/usr/bin/env python3

import math
import sys
import numpy
import scipy.special

# BORPStream failure probability calculator

# Parameters:
# N: the total number of real packets for the whole network
# r: (0 < r <= 1) the fraction of packets that are real packets;
#    the rest are dummies.  Therefore, there are N/r packets in total.
# f: the fanout factor at each MergeSplitNode (MSN), also called the
#    number of "flavours" each real packet entering an MSN can have;
#    each flavour of real packet corresponds to a particular one of the
#    f outputs of the MSN.  Dummy packets have no flavour.
# d: the depth of the MSN network.  There are d layers of f^{d-1} MSNs each.
#    The last layer outputs to f^d output buckets, each of size (N/r)/f^d.
# s: the internal size of the packet buffer at each MSN.

def flavour_dist(r, f, flavours_remaining=None, packets_remaining=None,
        prefix=None, scale=1.0):
    """Given the probability of a packet being real (r), and the total
    number of flavours (f), output a dictionary whose keys are the
    possible distributions of the number of real packets of each flavour
    (a tuple of length f) and whose values are the probability of that
    distribution of packets.  If flavours_remaining is not None, it is
    the number of remaining flavours to consider (so the keys of the
    output dictionary will have length flavours_remaining instead of f).
    If packets_remaining is not None, it is the number of packets
    remaining to be considered. If prefix is not None, it is prepended
    to each key.  If scale is not 1.0, it is multiplied by each
    probability value."""
    if prefix is None:
        prefix = ()
    if flavours_remaining is None:
        flavours_remaining = f
    if packets_remaining is None:
        packets_remaining = f
    if flavours_remaining == 0:
        return { prefix : scale * (1-r)**packets_remaining }
    retdict = {}
    for num_first_flavour in range(packets_remaining+1):
        retdict.update(flavour_dist(r, f, flavours_remaining - 1,
            packets_remaining - num_first_flavour,
            prefix + (num_first_flavour,),
            scale * scipy.special.comb(packets_remaining,
                    num_first_flavour, exact=True) *
                ((r/f) ** num_first_flavour)))
    return retdict

def transition_row(dist, capacity, curstate):
    """Given the output of flavour_dist, the capacity of the buffer, and
    a tuple representing the current state of the buffer (each element
    of the tuple is the number of packets of each flavour currently in
    the buffer), output a dictionary where the keys are the possible
    output states after processing f incoming packets and emitting one
    packet of each flavour (where possible).  The values are the
    probabilities of movine to that next state.  The special state
    "FAIL" means that the buffer capacity was exceeded."""
    retdict = {}
    for k in dist:
        if k == "FAIL" or curstate == "FAIL":
            # FAIL is an absorbing state; once you're in it, you're
            # stuck there
            newstate = "FAIL"
        else:
            newstate = tuple(map(lambda a,b: max(a+b-1,0), curstate+(0,), k))
            if sum(newstate) > capacity:
                newstate = "FAIL"
        if newstate != "FAIL":
            # Sort the entries
            newstate = tuple(sorted(newstate, reverse=True))
            # Remove the last (smallest) entry, since you can't increase
            # the capacity of the buffer unless at least one of the
            # flavours has no packets in the buffer.
            newstate = newstate[:len(curstate)]
        if newstate not in retdict:
            retdict[newstate] = dist[k]
        else:
            retdict[newstate] += dist[k]
    return retdict

def transition_matrix(f, r, capacity):
    """Given the number of flavours, the proportion of real packets, and
    the capacity of the buffer, produce the transition matrix for the
    action of receiving f packets and emitting one of each flavour
    (where possible).  Row/column 0 will be the start state (empty
    buffer), and row/column 1 will be the FAIL state where the buffer
    capacity has been exceeded.  The other row/column numbers are
    assigned to states arbitrarily."""

    statelist = [ (0,)*(f-1), "FAIL" ]
    statedict = { (0,)*(f-1): 0, "FAIL": 1 }

    M = []
    fdist = flavour_dist(r, f)
    statenum = 0
    # Note that statelist gets expanded as this loop executes, so we
    # definitely mean for len(statelist) to be checked every time
    # through the loop.
    while statenum < len(statelist):
        rowdict = transition_row(fdist, capacity, statelist[statenum])
        newrow = []
        for k in rowdict:
            # Has k been assigned a row/column id yet?
            if k not in statedict:
                kid = len(statelist)
                statelist.append(k)
                statedict[k] = kid
            else:
                kid = statedict[k]
            if len(newrow) <= kid:
                newrow += ([0.0] * (kid + 1 - len(newrow)))
            newrow[kid] += rowdict[k]
        M.append(newrow)
        statenum += 1

    # Now extend every row of M to the right length
    for i in range(statenum):
        M[i] += ([0.0] * (statenum - len(M[i])))

    return numpy.array(M)

def msn_failure_prob(num_packets, r, f, s):
    """Given the number of packets input to an MSN, the fraction r of
    them that are real, the number of flavours, and the size of the MSN
    internal buffer, return the probability that the buffer will
    overflow."""

    # We process the packets f at a time, so we need to reserve f-1
    # buffer slots to hold packets until we process them
    if s < f:
        return 1.0

    # The remaining capacity of the buffer
    capacity = s - (f-1)

    M = transition_matrix(f, r, capacity)

    # Each transition is for f packets, so we will invoke it
    # ceil(num_packets/f) times
    num_transitions = math.ceil(num_packets/f)

    Mpower = numpy.linalg.matrix_power(M,num_transitions)

    return Mpower[0][1]

def borpstream_failure_prob(N, r, f, d, s):
    """Given the total number of real packets (N), the proportion of all
    packets that are real (r), the fanout (f), depth (d), and buffer
    size (s), compute (an upper bound on) the overall network failure
    probability."""

    print("borpstream_failure_prob", N, r, f, d, s)
    sys.stdout.flush()
    # There will be N/r total packets, and each of the f^{d-1} MSNs in
    # each of d layers will receive (N/r)/f^{d-1} of them.
    msns_per_layer = f**(d-1)
    tot_msns = msns_per_layer * d

    per_msn_fail_prob = \
        msn_failure_prob(math.ceil((N/r)/msns_per_layer), r, f, s)

    # Union bound over the d*f^{d-1} MSNs
    tot_fail_prob = tot_msns * per_msn_fail_prob

    return tot_fail_prob

def borpstream_min_s(N, r, f, d, epsilon):
    """Given the total number of real packets (N), the proportion of all
    packets that are real (r), the fanout (f), depth (d), and the
    desired bound on the failure probability (epsilon), find the
    smallest s that satisfies the bound. Return s and its corresponding
    failure probability."""

    # s=1 is definitely too small
    s_lower = 1

    # Keep increasing s until we find one that works (but is probably not
    # minimal)
    while True:
        if s_lower < 20:
            s = s_lower + 8
        elif s_lower < 40:
            s = s_lower + 6
        else:
            s = s_lower + 2
        fp = borpstream_failure_prob(N, r, f, d, s)
        if fp <= epsilon:
            s_upper = s
            s_upper_fp = fp
            break
        else:
            s_lower = s

    # Now binary search to find the minimal s
    # Invariant: s_lower < s_min <= s_upper
    while s_upper > s_lower + 1:
        s = math.floor((s_upper+s_lower)/2)
        fp = borpstream_failure_prob(N, r, f, d, s)
        if fp <= epsilon:
            s_upper = s
            s_upper_fp = fp
        else:
            s_lower = s

    return s_upper, s_upper_fp

if __name__ == "__main__":
    epsilon = 2.0**(-80)
    if len(sys.argv) > 1:
        N = int(sys.argv[1])
    else:
        N = 1048576
    if len(sys.argv) > 2:
        r = float(sys.argv[2])
    else:
        r = 0.5
    if len(sys.argv) > 3:
        flist = [int(sys.argv[3])]
    else:
        flist = [2, 3, 4]
    if len(sys.argv) > 4:
        dlist = lambda f: [int(sys.argv[4])]
    else:
        dlist = lambda f: range(1, math.ceil(math.log(N/32)/math.log(f)))

    for f in flist:
        for d in dlist(f):
            if f**d > N/32:
                continue
            s, fp = borpstream_min_s(N, r, f, d, epsilon)
            if fp > 0.0:
                lfp = math.log(fp)/math.log(2.0)
            else:
                lfp = "-inf"
            print(N, f, d, s, fp, lfp)
