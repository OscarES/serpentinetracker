class epics :
    def __init__(self, so) :
        # keep local reference of serpentine object
        self.so = so

        # list of pv and variable associations
        self.ass = []        

    def buildIOC(self, appName = 'test', iocName='instance') :
        # loop over elements and create cmd files        
        for e in self.so.beamline :
            print e.name, e.__class__.__name__ 

    def startIOC(self, so) :
        pass 

    def writeIOC(self) :
        pass

    def readIOC(self) :
        pass


