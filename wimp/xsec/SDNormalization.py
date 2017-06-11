import numpy as np

class SDNormalization:
    
    def __init__(self):
        self.mp = 938.272 * units.MeV
        self.mn = 939.565 * units.MeV
        self.mup = 5.586 # in nuclear magnetons
        self.mun = -3.826 
        self.gp = 1
        self.gn = 1

    def single_particle(self,nucleus,Mx,to_proton=True):
        m = nucleus.mass()
        mu = m * Mx / (Mx + m)
        muN = Mx * self.mp / (self.mp + Mx)
        if not to_proton:
            muN = Mx * self.mn / (self.mn + Mx)
 
        # Odd type:
        g2 = self.gp*self.gp
        if n.N()%2==1:
          g2 = self.gn*self.gn

        g2 *= ( nucleus.J()*(nucleus.J()+1) + 0.75 -\
                nucleus.L()*(nucleus.L()+1))**2 /\
               (4*nucleus.J()*(nucleus.J()+1) )

        g2N = self.gp**2 * 0.75 #J(J+1)

        if not to_proton:
            muN = Mx * self.mn / (self.mn + Mx)
            g2N = self.gn**2 * 0.75
    
        return muN*muN * g2N / ( mu*mu*g2)
        
    def odd_group(self,nucleus,Mx,to_proton=True):
        m = nucleus.mass()
        mu = m * Mx / (Mx + m)
        muN = Mx * self.mp / (self.mp + Mx)
        if not to_proton:
            muN = Mx * self.mn / (self.mn + Mx)

        # Odd type:
        g2 = self.gp*self.gp * (nucleus.J()+1)/nucleus.J()
        g2 *= ((nucleus.magnetic_moment() - nucleus.J())/(self.mup - 1))**2
        if n.N()%2==1:
            g2 = self.gn*self.gn * (nucleus.J()+1)/nucleus.J()
            g2 *= (nucleus.magnetic_moment()/self.mun)**2
 
        g2N = self.gp**2 * 0.75 #J(J+1)

        if not to_proton:
            muN = Mx * self.mn / (self.mn + Mx)
            g2N = self.gn**2 * 0.75
    
        return muN*muN * g2N / ( mu*mu*g2)
        
    def spin_content(self,nucleus,Mx,to_proton=True):
        m = nucleus.mass()
        mu = m * Mx / (Mx + m)
        muN = Mx * self.mp / (self.mp + Mx)
        if not to_proton:
            muN = Mx * self.mn / (self.mn + Mx)


        g2 = (nucleus.J()+1) / nucleus.J()
        g2 *= (self.gp * nucleus.proton_spin() + \
               self.gn*nucleus.neutron_spin())**2

        g2N = self.gp**2 * 0.75 #J(J+1)

        if not to_proton:
            muN = Mx * self.mn / (self.mn + Mx)
            g2N = self.gn**2 * 0.75
    
        return muN*muN * g2N / ( mu*mu*g2)
        

    def normalization(self,nucleus,Mx,method='spin_content',to_proton=True):
        if method== 'odd_group':
            return self.odd_group(nucleus,Mx,to_proton)
        elif method == 'single_particle':
            return self.single_particle(nucleus,Mx,to_proton)

        return self.spin_content(nucleus,Mx,to_proton)
