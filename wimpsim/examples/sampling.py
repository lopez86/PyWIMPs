from ROOT import TH1F
from ROOT import TCanvas
from ROOT import kRed, kBlue, kGreen
from ROOT import TF1

from ..astro import AstroModel
from ..xsec import InteractionModel
from ..sampling import UniformWeightedSampler
from ..sampling import MaxwellWeightedSampler
from ..sampling import AcceptRejectSampler
from ..sampling import MCMCSampler
from .. import units

import numpy as np

def CompareMethods():
    h1 = TH1F('unif','Uniform;E_{r} [keV];Rate [shape only]',100,0,200)
    h2 = TH1F('maxw','Maxwell;E_{r} [keV];Rate [shape only]',100,0,200)
    h3 = TH1F('ar','A/R;E_{r} [keV];Rate [shape only]',100,0,200)
   



    h1.Sumw2()
    h2.Sumw2()
    h3.Sumw2()

    am = AstroModel()
    im = InteractionModel()
    unif = UniformWeightedSampler(am,im)
    maxw = MaxwellWeightedSampler(am,im)
    ar = AcceptRejectSampler(am,im)

    unif.initialize()
    maxw.initialize()
    ar.initialize()

    C = units.speed_of_light * np.sqrt(0.5 / im.cross_section.Mt * units.keV)  * (im.cross_section.Mt+im.cross_section.Mx) / im.cross_section.Mx
    f = TF1('myfun','[0] * (TMath::Erf( ([1] * sqrt(x) + [2]) / [3] ) - TMath::Erf( ([1]*sqrt(x)-[2])/[3]))',0,200)
    f.SetNpx(1000)
    f.FixParameter(1,C)
    f.FixParameter(2,np.sqrt(am.velocity.vE.dot(am.velocity.vE)))
    f.FixParameter(3,am.velocity.v0)
    f.SetParameter(0,1)   

    for i in range(100000):
        if i%10000==0:
          print('Getting event: '+str(i))
        us = unif.sample()
        ms = maxw.sample()
        ars = ar.sample()

        h1.Fill(us.Er / units.keV,us.weight)
        h2.Fill(ms.Er / units.keV,ms.weight)
        h3.Fill(ars.Er / units.keV,ars.weight)

    h1.Scale(1./h1.GetBinContent(1))
    h2.Scale(1./h2.GetBinContent(1))
    h3.Scale(1./h3.GetBinContent(1))

    h1.SetLineColor(kRed)
    h2.SetLineColor(kBlue)
    h3.SetLineColor(kGreen+1)
    h1.SetMarkerColor(kRed)
    h2.SetMarkerColor(kBlue)
    h3.SetMarkerColor(kGreen+1)

    a = f.Eval(0)
    f.SetParameter(0, 1. / a)

    f.SetLineWidth(2)
    h1.SetLineWidth(2)
    h2.SetLineWidth(2)
    h3.SetLineWidth(2)

    return h1,h2,h3,f

def main():
    CompareMethods()

if __name__=='__main__':
    main()


