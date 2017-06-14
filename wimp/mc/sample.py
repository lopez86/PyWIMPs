""" sample.py 

    Module holding a class representing a single event.

"""
__author__ = "Jeremy P. Lopez"
__date__      = "June 2017"
__copyright__ = "(c) 2017, Jeremy P. Lopez"

class Sample:
    """ A class holding the information from a single
        event ('Sample') in the detector.
 
        Attributes:
            Er_reco: Reconstructed recoil energy
            recoil_vec_reco: Recoil direction vector
            gen_weight: Generator weight
            det_weight: Detector weight (default: 1)
            wimp_vec: Initial WIMP velocity vector
    """
    def __init__(self,Er,vrecoil,weight,vwimp):
        """ Initialize. Reconstructed variables 
            are initialized to the true values 
            with a detector weight of 1.

            Args:
                Er: Recoil energy
                vrecoil: Recoil direction vector
                weight: Generator weight
                vwimp: WIMP velocity vector
        """

        self.Er = Er
        self.recoil_vec = vrecoil
        self.gen_weight = weight
        self.wimp_vec = vwimp
        self.det_weight = 1.0

    @property
    def weight(self):
        """ The total weight (product of generator and
            detector level weights)
        """
        return self.gen_weight * self.det_weight
 
    @property
    def Er(self):
        """ The recoil energy. """
        return self._Er     

    @Er.setter
    def set_Er(self,er):
        """ Set the recoil energy.

            Args:
                er: The recoil energy
        """
        self._Er = er
        self.Er_reco = er

    @property
    def recoil_vec(self):
        """ The recoil direction vector. """
        return self._recoil_vec

    @recoil_vec.setter
    def set_recoil_vec(self,rv):
        """ Set the recoil direction vector.

            Args:
                rv: The recoil direction
        """
        self._recoil_vec = rv
        self.recoil_vec_reco = rv
