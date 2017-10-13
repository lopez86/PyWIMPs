""" maximum_gap.py

Maximum gap method described by Yellin. This is a
method for finding an upper limit with an unknown background
and a non-zero number of events.

TODO:
    Complete module

"""



def maximum_gap_ul(E,Emin,Emax,model,cl=0.9,Ntrials=1000000):
    print('This is still under development! Returning')
    return -1

    Evec = set(E + [Emin,Emax])
    gaps = np.zeros(len(Evec)-1)

    model.initialize()
    for i in range(Ntrials):
        pass             


