from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import collections
import numpy as np
import scipy.integrate as integrate
from . import dusts
__ver__ = '1.0'

class sfh_wrapper(object):
    """ sfh_wrapper class.  EzGal wraps this class around the sfh function.
    It takes care of the details of passing or not passing parameters """

    func = ''       # sfh function
    args = ()       # extra arguments to pass on call
    has_args = False    # whether or not there are actually any extra arguments

    def __init__(self, function, args):
        """ wrapper_obj = ezgal.sfhs.wrapper(function, args)

        wrapper class.  EzGal wraps this class around the sfh function.
        It takes care of the details of passing or not passing parameters """

        self.func = function

        if type(args) == tuple and len(args) > 0:
            self.has_args = True
            self.args = args

    def __call__(self, val):

        if self.has_args:
            return self.func(val, *self.args)
        else:
            return self.func(val)

class numeric(object):
    ages = np.array([])
    sfr = np.array([])

    def __init__(self, ages, sfr):
        """ numeric_obj = ezgal.sfhs.numeric(ages, sfrs)

        wrapper class for making a numeric star formation history callable.
        Pass a list of ages and relative star formation rates.  Ages should be in gyrs. """

        self.ages = np.asarray(ages)
        self.sfr = np.asarray(sfr)

    def __call__(self, val):
        return np.interp(val, self.ages, self.sfr)

def exponential(t, tau):
    """ ezgal.sfhs.exponential(ages, tau)

    exponentially decaying star formation history with
    e-folding time scale of tau gyrs """

    return np.exp(-1.0*t/tau)

def constant(t, length):
    """ ezgal.sfhs.constant(ages, length)

    Burst of constant starformation from t=0 to t=length """

    if type(t) == np.ndarray:
        sfr = np.zeros(t.size)
        m = t <= length
        if m.sum(): sfr[m] = 1.0
        return sfr
    else:
        return 0.0 if t > length else 1.0

def exponential_truncation(t, tau, t_cut, tau_cut):
    """ ezgal.sfhs.exponential_truncation(ages, tau, age_cut, tau_cut)

    Exponentially decaying star formation with truncation """

    sfr = np.exp(-1.0*t/tau)
    if type(t) == np.ndarray:
        index = t > t_cut
        sfr[index] = np.exp((t_cut -t[index])/tau_cut) * np.exp(-1.0*t_cut/tau)
    else:
        if t > t_cut:
            sfr = sfr * np.exp((t_cut - t)/tau_cut)
    return sfr

