import numpy as np
import pickle
import pandas as pd

def pyko_to_normal_panda(pkodata):
    df = pd.DataFrame({
            "j"    : pkodata.j.magnitude,
            "stepn" : pkodata.stepn.magnitude,
            "time" : pkodata.time.magnitude,
            "mat" : pkodata.mat.magnitude,
            "pos" : pkodata.pos.magnitude,
            "dr" : pkodata.dr.magnitude,
            "rho0" : pkodata.rho0.magnitude,
            "rho" : pkodata.rho.magnitude,
            "up" : pkodata.up.magnitude,
            "ie" : pkodata.ie.magnitude,
            "pres" : pkodata.pres.magnitude,
            "mass" : pkodata.mass.magnitude*4*np.pi, # correction needed for spherical true mass
            "temp" : pkodata.temp.magnitude,
            "cs" : pkodata.alocal.magnitude,
            "phase" : pkodata.phase.magnitude,
            "etot" : pkodata.etot.magnitude,
            "dtminj" : pkodata.dtminj.magnitude,
            "ent" : pkodata.entropy.magnitude,
            })
    return df

class plume_class:
    """Class to hold plume scaling data."""  # this is a documentation string for this class
    def __init__(self): # self is the default name of the object for internal referencing of the variables in the class
        """A function to initialize the class object.""" # this is a documentation string for this function
        self.rplumeinitarr = []   
        self.pinitarr    = []
        self.velinitarr = []
        self.einitarr  = []
        self.minitarr = []
        self.rhoinitarr = [] 
        self.rstallarr = []
        self.tstallarr = []
        self.vkearr = []
        self.labelneb  = []
        self.pneb  = []
        self.rhonebarr = []
        self.gammanebarr = []
        self.csnebarr = []
        self.frhonebref = []
        self.symarr = []
        self.cauchy = []
        self.pimass = []
        self.pitime = []
        self.piradius = []
#
plume = plume_class()
#
# varying simulation conditions
plume.rplumeinitarr = np.asarray([1.e3,5.e3,25.e3,125.e3]) # m
plume.pinitarr = np.asarray([100.e9,150.e9,200.e9,250.e9,350.e9,450.e9]) #Pa
plume.velinitarr = np.asarray([0.]) # m/s

plume.labelneb = ['0.01 Pa','0.1 Pa','1 Pa']
plume.pneb = ['p01Pa','p1Pa','1Pa']
plume.rhonebarr = np.asarray([1.406419e-11,1.406419e-10,1.406419e-09])*1000. # kg/m3
plume.presnebarr = np.asarray([0.01,0.1,1])
plume.gammanebarr = np.asarray([1.4,1.4,1.4]) # dimless
plume.csnebarr = np.sqrt(plume.gammanebarr*plume.presnebarr/plume.rhonebarr) # m/s
plume.rhonebref = 1.406419e-11*1000. # kg/m3

plume.symarr = ['.','.','.']

nrp = len(plume.rplumeinitarr) # varying plume size
npp = len(plume.pinitarr) # varying plume pressure
nvp = len(plume.velinitarr) # varying initial velocity
npn = len(plume.pneb) # varying nebula properties; nebula pressure or nebular dust

plume.einitarr = np.zeros([npp,nrp])
plume.minitarr = np.zeros([npp,nrp])
plume.rhoinitarr = np.zeros([npp,nrp])
plume.rstallarr = np.zeros([npp,nrp,npn])
plume.tstallarr = np.zeros([npp,nrp,npn])
plume.cauchy = np.zeros([npp,nrp,npn])
plume.pimass = np.zeros([npp,nrp,npn])
plume.pitime = np.zeros([npp,nrp,npn])
plume.piradius = np.zeros([npp,nrp,npn])
plume.vkearr = np.zeros([npp,nrp])

for ineb in range(npn): # vary nebula
    fout='./data/silica/vp-silica-grid-'+plume.pneb[ineb]+'/vp-silica-sphere-grid-'
    for ipp in range(npp): # vary plume pressure
        for irp in range(nrp): # vary plume radius
            for ivel in range(nvp):       # vary plume velocity     
                fileid = 'P'+str(np.round(plume.pinitarr[ipp]/1.e9))+'-R'+str(np.round(plume.rplumeinitarr[irp]/1.e3))+'-V'+str(np.round(plume.velinitarr[ivel]/1.e3))
                outputfilename = fout+fileid+'.dat'
                stalloutputfilename = fout+fileid+'.dat-stall.one'
                print('outputfile',outputfilename)
                with open(outputfilename,"rb") as f:
                    pkodata = pickle.load(f) # keeps units
                with open(stalloutputfilename,"rb") as f:
                    stallpkodata = pickle.load(f) # keeps units
                if 0:
                    # print units
                    print('pyKO output file units are the same as the input file units:')
                    print('   Time        ',pkodata.time.units)
                    print('   Position    ',pkodata.pos.units)
                    print('   Density     ',pkodata.rho.units)
                    print('   Part. vel.  ',pkodata.up.units)
                    print('   Int. energy ',pkodata.ie.units)
                    print('   Mass        ',pkodata.mass.units)
                    print('   Temperature ',pkodata.temp.units)
                    print('   Sound speed ',pkodata.alocal.units)
                    print('   Pressure    ',pkodata.pres.units)
                    print('   Stress      ',pkodata.sigmar.units)
                    print('   Sp. Entropy ',pkodata.entropy.units)
                #jmaxmat1=max(np.where(pkodata.mat.magnitude==1)[0])
                pko = pyko_to_normal_panda(pkodata)
                stallpko = pyko_to_normal_panda(stallpkodata)
                imat1 = np.where((pko['mat'] == 1))[0]
                etotal = np.sum(pko['ie'][imat1]*pko['mass'][imat1])
                mtotal = np.sum(pko['mass'][imat1])
                #print('Etotal (J) = ',etotal)
                plume.einitarr[ipp,irp]=etotal
                plume.minitarr[ipp,irp]=mtotal
                plume.rhoinitarr[ipp,irp]=pko['rho'][imat1[0]]
                plume.vkearr[ipp,irp] = np.sqrt(etotal*2/mtotal)
                # stall info
                plume.rstallarr[ipp,irp,ineb] = np.max(stallpko['pos'][np.max(imat1)])
                plume.tstallarr[ipp,irp,ineb] = stallpko['time'][0]
                # cauchy number
                plume.cauchy[ipp,irp,ineb] = (etotal*2/mtotal)/(plume.csnebarr[ineb]*plume.csnebarr[ineb])
                mstall = (4/3*np.pi)*plume.rhonebarr[ineb]*np.power(plume.rstallarr[ipp,irp,ineb],3)
                plume.pimass[ipp,irp,ineb] = mstall/mtotal
                plume.pitime[ipp,irp,ineb] = plume.vkearr[ipp,irp]*plume.tstallarr[ipp,irp,ineb]/plume.rplumeinitarr[irp]
                plume.piradius[ipp,irp,ineb] = plume.rstallarr[ipp,irp,ineb] / plume.rplumeinitarr[irp]

