from elements import Twiss, Jitter

def laodConfig(filename=''):

    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    t = Twiss()
    j = Jitter()

    for (i,l) in enumerate(lines):
        if l.startswith('#') or l.isspace() or not l: continue
        ls = l.split()
        if ls[0] in t.__dict__.keys():
            t.__setattr__(ls[0],float(ls[-1]))
        elif ls[0] in j.__dict__.keys():
            j.__setattr__(ls[0],float(ls[-1]))
        else:
            print 'TwissLoader.loadTwiss > line '+str(i)+' ignored: no recognised parameter'

    return t,j      
    
