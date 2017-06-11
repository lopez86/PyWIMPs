class SINormalization:
    def __init__(self):
        self.mp = 938.272 * units.MeV
        self.mn = 939.565 * units.MeV
        self.gp = 1
        self.gn = 1

    def normalize(self,nucleus,Mx,o_proton=True):
        m = nucleus.mass()
        mu = m * Mx / (Mx + m)
        muN = Mx * self.mp / (self.mp + Mx)
        if not to_proton:
            muN = Mx * self.mn / (self.mn + Mx)

        g = nucleus.Z() * self.gp + nucleus.N() * self.gn
        gN = self.gp

        if not to_proton:
            muN = Mx * self.mn / (self.mn + Mx)
            gN = self.gn
    
        return muN*muN * gN * gN / ( mu*mu*g*g)
        
