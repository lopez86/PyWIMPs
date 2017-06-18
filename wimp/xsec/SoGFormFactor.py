""" SoGFormFactor.py

Sum of Gaussians form factor. This is the form factor from a 
charge density defined by the sum of a series of Gaussian
distributions with different means but identical scales.

Data defined for:
C-12, O-16, Si-28, S-32, Ca-40, Pb-208
arXiv:0608035 (hep-ph)
"""

__author__ = 'Jeremy P. Lopez'
__date__ = 'June 2017'
__copyright__ = '(c) 2017, Jeremy P. Lopez'


import numpy as np
from .. import units
from .FormFactor import FormFactor
_u = units.fm / units.hbarc

SoGRi = {
         (6,12):np.array((0.0,0.4,1.0,1.3,1.7,2.3,2.7,3.5,4.3,5.4,6.7)) \
                * _u,
         (8,16):np.array((0.4,1.1,1.9,2.2,2.7,3.3,4.1,4.6,5.3,5.6,5.9,6.4)) \
                * _u,
         (14,28):np.array((0.4,1.0,1.9,2.4,3.2,3.6,4.1,4.6,5.1,5.5,6.0,6.9)) \
                 * _u,
         (16,32):np.array((0.4,1.1,1.7,2.5,3.2,4.0,4.6,5.0,5.5,6.3,7.3,7.7)) \
                 * _u,
         (20,40):np.array((0.4,1.2,1.8,2.7,3.2,3.6,4.3,4.6,5.4,6.3,6.6,8.1)) \
                 * _u,
         (82,208):np.array((0.1,0.7,1.6,2.1,2.7,3.5,4.2,5.1,6.0,6.6,7.6,8.7)) \
                 * _u
        }

SoGQi = {
         (6,12):np.array((0.016690,0.050325,0.138621,0.180515,0.219097,0.278416,
                 0.058779,0.057817,0.007739,0.02001,0.00007)),
         (8,16):np.array((0.057056,0.195701,0.311188,0.224321,0.059946,0.135714,
                 0.000024,0.013961,0.000007,0.000002,0.002096,0.000002)),
         (14,28):np.array((0.033149,0.106452,0.206866,0.286391,0.250448,0.056944,
                  0.016829,0.039630,0.000002,0.000938,0.000002,0.002366)),
         (16,32):np.array((0.045356,0.067478,0.172560,0.324870,0.254889,0.101799,
                  0.022166,0.002081,0.005616,0.000020,0.000020,0.003219)),
         (20,40):np.array((0.042870,0.056020,0.167853,0.317962,0.155450,0.161897,
                  0.053763,0.032612,0.004803,0.004541,0.000015,0.002218)),
         (82,208):np.array((0.003845,0.009724,0.033093,0.000120,0.083107,0.080869,
                   0.139957,0.260892,0.336013,0.033637,0.018729,0.000020))
        }


SoGData = {
         (6,12):(2.469*_u,SoGRi[(6,12)],SoGQi[(6,12)],1.20*_u),
         (8,16):(2.711*_u,SoGRi[(8,16)],SoGQi[(8,16)],1.30*_u),
         (14,28):(3.121*_u,SoGRi[(14,28)],SoGQi[(14,28)],1.30*_u),
         (16,32):(3.258*_u,SoGRi[(16,32)],SoGQi[(16,32)],1.35*_u),
         (20,40):(3.480*_u,SoGRi[(20,40)],SoGQi[(20,40)],1.45*_u),
         (82,208):(5.503*_u,SoGRi[(82,208)],SoGQi[(82,208)],1.70*_u),
        }


class SoGFormFactor(FormFactor):
    """ The sum of Gaussians form factor.

        This comes from a charge distribution defined by a sum of
        p(r,k) = a(k)[exp(- (r - R(k))^2/g^2] + exp(-(r+R(k))^2/g^2 
        over an integer k, with p(r) = 0 
        
        where a(k) are factors related to some fractional charge Q(k)
        and the radius scale R(k)

        The data should be packaged in such a way that
        data[1] = numpy array of R(k) 
        data[2] = numpy array of Q(k)
        data[3] = g * sqrt(1.5)
        R(k) and g should be in units of inverse energy.

    """

    def __init__(self,data=SoGData[(6,12)]):
        """ Initialize with carbon as the default. 
 
            Args:
                data: The nuclear data needed for this calculation
                      Default: Carbon-12
                      Given in units of inverse energy (length/ (hbar*c) )

        """
        self._data = data
        self.F0 = 1
        self.random = np.random

    def initialize(self):
        """ Get the normalization constant. """
        self.F0 = 1
        f0 = self.ff2(0)
        self.F0 = f0
 

    def set_data(self,data):
        """ Set the data.

            Args:
                data: The nuclear data needed for this calculation
                      Given in units of inverse energy (length/ (hbar*c) )

        """
        self._data = data

    def set_params(self,pars):
        """ Set the parameters from a dictionray

            Args:
                pars: {string}
 
            Parameters:
                SoGData: The data in the correct format and in units
                        of inverse energy.
        """

        if 'SoGData' in pars.keys():
            self._data = pars['SoGData']

    def ff2(self,Q2):
        """ Calculate the form factor as a function of Q^2. 

            Args:
                Q2: The squared momentum transfer (a positive value)
        """
        gamma = _data[3] / np.sqrt(1.5) 
        A = exp(-0.25*Q2*gamma*gamma)
        q = np.sqrt(Q2)

        F = 0
        for i in range(1,len(_data[1].size+1)):
            F += _data[2][i-1] / (1 + 2 *(_data[1][i-1]/gamma)**2)  \
                 * (np.cos(q * self._data[1][i-1]) 
                    + 2 * (self._data[1][i-1]/gamma)**2
                    * np.sin(q*self._data[1][i-1]) / (q * self._data[1][i-1]))

        return F*F /self.F0 / self.F0
