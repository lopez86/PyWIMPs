class Nucleus:
    def __init__(self,pars=None):
        if pars is not None:
            self.set_params(pars)
        else:
            self._A = 0
            self._Z = 0
            self._mass = 0
            self._mag_mom = 0
            self._spin = 0
            self._Sp = 0
            self._Sn = 0

    def A(self):
        return self._A
    def Z(self):
        return self._Z
    def N(self):
        return self._A - self._Z
    def mass(self):
        return self._mass
    def magnetic_moment(self):
        return self._mag_mom
    def spin(self):
        return self._spin
    def proton_spin(self):
        return self._Sp
    def neutron_spin(self):
        return self._Sn
