class Config:
    def __init__(self, fname=None):
        import ConfigParser
        
       # Default values, see explanations below:        
        configDic = {
        'isMetallic': 'True',
        'numProc': '1',
        'paraPrefix': '',
        'paraPostfix': '',
        'pwscfInput': 'scf.in',
        'pwscfOutput': 'scf.out',
        'phInput': 'ph.in',
        'phOutput': 'ph.out',
        'dynmatInput': 'dynmat.in',
        'dynmatOutput': 'dynmat.out',
        'q2rInput': 'q2r.in',
        'q2rOutput': 'q2r.out',
        'matdynInput': 'matdyn.in',
        'matdynOutput': 'matdyn.out',
        'matdynModes': 'matdyn.modes'
        }

        try:
          if fname == None:
             raise NameError("Config should be initialized with a filename")
        except NameError:
            raise
        
        self.config = ConfigParser.SafeConfigParser(configDic)
        self.config.read(fname)
        
        # all the relevant input files must be preconfiguered for specific tasks 
        # before using this class

        # parallelization parameters
        self.numProc = self.config.getint('Config','numProc')
        self.paraPrefix = self.config.get('Config', 'paraPrefix')
        self.paraPostfix = self.config.get('Config', 'paraPostfix')  
        
        self.pwscfInput = self.config.get('Config', 'pwscfInput')
        # pwscf output file relevant to 'total energy' as well as 'geometry' tasks           
        self.pwscfOutput = self.config.get('Config', 'pwscfOutput')
        
        self.phInput = self.config.get('Config', 'phInput')
        self.phOutput = self.config.get('Config', 'phOutput')
        
        # dynmat input/output file relevant to 'single phonon' task        
        self.dynmatInput = self.config.get('Config', 'dynmatInput')
        self.dynmatOutput = self.config.get('Config', 'dynmatOutput')
        
        # input/output files relevant to 'multiple phonon' task    
        self.q2rInput = self.config.get('Config', 'q2rInput')
        self.q2rOutput = self.config.get('Config', 'q2rOutput')            
        self.matdynInput = self.config.get('Config', 'matdynInput')
        self.matdynOutput = self.config.get('Config', 'matdynOutput')
        self.matdynModes = self.config.get('Config', 'matdynModes')     
        