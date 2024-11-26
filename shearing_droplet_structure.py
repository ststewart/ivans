#SJL 3/18
#script to calculate the equilibrium size of droplets swept up in a bow-shock

###########################################################
###########################################################
###########################################################
import numpy as np
import sys
import math

import struct

from scipy import constants as const
from scipy.optimize import fsolve
from scipy.optimize import newton
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

###########################################################
###########################################################
###########################################################
#CONST


###########################################################
###########################################################
###########################################################
#Define a structure which contains all the information for a shearing droplet

class shearing_droplet(object):
    ###################################################
    ###################################################
    #takes the mass of the two bodies as inputs
    def __init__(self):

        #given droplet properties
        self.flag_cond=0 #whether the condensate is SiO2 (0) or Fe (1)
        self.Rcond=0.0 #radius
        self.rho_cond=0.0 #density
        self.Tcond=0.0 #temperature
        self.cs_cond=0.0 #Sound speed in condensate
        self.comp_cond={}
        self.flag_comp_mol=0
        self.flag_visc_cond_model=0

        #properties for aggregates
        self.Vfrac_ag=0.5
        self.rho_mono_ag=3E3
        self.r_mono_ag=1E-6

        #derived droplet parameters
        self.gamma=0.0 #surface tension
        self.m=0.0 #mass of condensate
        self.visc_cond=0.0

        #given vapor properties
        self.flag_vap=0 #flag as to whether water is H/He (0) or H2O (1)
        self.vap_prop=[] #array of properties needed for vapor
        self.rho_vap=0.0
        self.Tvap=0.0
        self.vvap=0.0
        self.cs=np.nan #if still nan then derived as ideal gas

        #derived material properties
        self.visc=0.0
        self.sigma_coll=0.0
        self.cpcv=0.0
        self.ma=0.0
        self.pvap=0.0
        self.MFP=0.0

        #dynamic parameters
        self.prefac=4.0 #prefactor for calculating balance between drag and 
        self.flag_CD_form_dynamic=1 #flag as to what CD form to use if in dynamic regime
        self.flag_CD_self_consistent=1 #flag as to whether to force the CD_form to be self-consistent (e.g., if in Epstein, use epstein)

        #derived dynamic properties
        self.vcrit=0.0 #critical velocity for breakup
        self.Req=0.0 #equilibrium size of a droplet with these gas properties based on force balance
        self.flag_FB_stable=True #whether the droplet is stable by force balance
        self.Rbreakup_Wec=np.asarray([]) #radius corresponding to the critical we for RT instability
        self.Rbreakup_Wec_inv=np.asarray([]) #radius corresponding to the critical we for RT instability assuming zero viscosity
        self.Rbreakup_Wec2=np.asarray([]) #radius corresponding to the critical we for RT instability
        self.tbreakup_cs=0.0 #droplet breakup time based on sound speed 
        self.tbreakup=0.0 #droplet breakup time based on experimental scalings
        self.CD=2.0
        self.Fdrag=0.0

        #Flag as to which regime the system is in 0: Epstein, 1: Dynamic
        self.flag_drag_regime=1

        #which regime of breakup are you in
        self.flag_breakup_regime=0
        self.regimeII_prefac=10.
        
        #Non-dimensional numbers
        self.Re=0.0 #Reynolds number
        self.Wec=0.0 #critical weber number for RT breakup for droplet
        self.Wec2=0.0 #critical weber number for SEI breakup for droplet

        #Normalizing times, lengthscales etc.
        tv_plus=0.0
        tI_plus=0.0

        #stopping time properties
        self.toll_stop=0.0
        self.Ncalc_stop=0
        self.tstop=0.0
        self.Lstop=0.0
        self.tcalc_stop=np.empty(0)
        self.sol_stop=np.empty([1,4])
        self.flag_stop_stab=1
        
        #properties if considering a second wave
        self.flag_hit_2nd_wave=0
        self.tstop1=0.0
        self.Lstop1=0.0
        self.tstop2=0.0
        self.Lstop2=0.0
        self.deltav_L=None
        self.flag_vap2=None
        self.rho_vap2=None
        self.Tvap2=None
        self.deltav2=None

    ###################################################
    ###################################################
    #function to intialize droplet
    def init_droplet(self,flag_cond,Rcond,Tcond,flag_vap,vap_prop,rho_vap,Tvap,vvap,cs,prefac,flag_CD_form_dynamic,flag_CD_self_consistent=1,cs_cond=1E3,regimeII_prefac=10.,comp={"SiO2": 100},flag_comp_mol=0,flag_visc_cond_model=0,Vfrac_ag=0.5, rho_mono_ag=3E3, r_mono_ag=1E-6,deltav_L=None,flag_vap2=None,rho_vap2=None,Tvap2=None,deltav2=None):

        ########
        #record the given droplet information
        ########
        self.flag_cond=flag_cond
        self.Rcond=Rcond
        self.Tcond=Tcond
        self.cs_cond=cs_cond
        self.comp=comp
        self.flag_visc_cond_model=flag_visc_cond_model

        if self.flag_cond==2: #dust bunny (i.e., agregate)
            self.Vfrac_ag=Vfrac_ag
            self.rho_mono_ag=rho_mono_ag
            self.r_mono_ag=r_mono_ag

        
        ########
        #record the shear properties
        ########
        self.prefac=prefac

        ########
        #derived condensate properties
        ########
        self.calc_gamma()
            
        if self.flag_cond==0: #silicate
            self.rho_cond=3E3
        elif self.flag_cond==1: #metal
            self.rho_cond=7.0E3
        elif self.flag_cond==2: #dust bunny (i.e., agregate)
            self.rho_cond=self.Vfrac_ag*self.rho_mono_ag
       
        else:
            raise ValueError('Unknown flag_cond')
        
        self.calc_visc_cond()

        self.m=4.0*np.pi*self.rho_cond*(self.Rcond**3)/3.0

        ########
        #record the vapor properties
        ########
        self.flag_vap=flag_vap
        self.vap_prop=vap_prop
        self.rho_vap=rho_vap
        self.Tvap=Tvap
        self.vvap=vvap
        self.cs=cs

        ########
        #derived vapor properties
        ########
        #cpcv and ma
        if self.flag_vap==0:
            self.cpcv=1.41
            XH2=self.vap_prop[0]
            self.ma=(2*XH2+4*(1-XH2))*1E-3
            self.sigma_coll=2.0E-19 #Tabata + Shirai 2000
        elif self.flag_vap==1:
            self.cpcv=1.33
            self.ma=18E-3
            self.sigma_coll=140E-20 #Mammali et al. 2016
        else:
            raise ValueError('Unknown flag_vap')

        #if no sound speed given calculate as an ideal gas
        if np.isnan(cs):
            cs=np.sqrt(self.cpcv*const.R*self.Tvap/(self.ma))
        self.cs=cs

        #calculate the viscosity
        self.calc_visc()

        #pressure from ideal gas
        self.pvap=self.rho_vap*const.R*self.Tvap/self.ma

        #mean free path
        self.MFP=const.k*self.Tvap/(np.sqrt(2)*self.sigma_coll*self.pvap)
        
        #########
        #2-step evolution properties
        self.deltav_L=deltav_L
        self.flag_vap2=flag_vap2
        self.rho_vap2=rho_vap2
        self.Tvap2=Tvap2
        self.deltav2=deltav2

        ########
        #check drag regime
        ########
        self.flag_CD_form_dynamic=flag_CD_form_dynamic
        self.flag_CD_self_consistent=flag_CD_self_consistent
        
        self.calc_drag_regime()
        
        if self.flag_CD_self_consistent==1:
            if self.flag_drag_regime==0:
                self.flag_CD_form=0
            else:
                self.flag_CD_form=self.flag_CD_form_dynamic
        else:
            self.flag_CD_form=self.flag_CD_form_dynamic


        ########
        #calc ND numbers etc.
        ########
        self.We=self.rho_vap*self.vvap**2*2*self.Rcond/self.gamma
        self.Re=self.calc_Re(self.Rcond,self.vvap)
        self.Oh=self.visc_cond/np.sqrt(self.rho_cond*self.gamma*2*self.Rcond)

        ########
        #calc drag 
        ########
        self.CD=self.calc_CD(self.Rcond,self.vvap,flag_CD_form=self.flag_CD_form)
        self.Fdrag=0.5*self.CD*np.pi*(self.Rcond**2)*self.rho_vap*(self.vvap**2)

        ########
        #calc timescales
        ########
        self.tv_plus=(self.visc_cond/self.rho_cond/((self.Fdrag/self.m)**2))**(1./3)
        self.tI_plus=1.0*(2.0*self.Rcond/self.vvap)*np.sqrt(self.rho_cond/self.rho_vap)

        ########
        #calculate the critical values
        ########
        #based on force balance
        self.calc_vcrit()
        self.Req=self.calc_Req(self.vvap,flag_CD_form=self.flag_CD_form)[0]

        #based on experimentally constrained scalings
        self.regimeII_prefac=regimeII_prefac
        self.calc_breakup_Wec(self.vvap,self.Rcond)

        

    ###################################################
    ###################################################
    #function to calculate the surface tension of the condensate
    #data from 
    def calc_gamma(self):

        #SiO2: data from Boca 2003
        if self.flag_cond==0:
            a=251.7E-3
            b=-0.0310E-3
            self.gamma=a-b*(self.Tcond-273.15)
            #self.gamma=0.5*self.gamma

        #Fe: data from
        elif self.flag_cond==1:
            TL=1538.
            gammaL=1.92
            gammaT=-3.97E-4
            self.gamma=gammaL+gammaT*(self.Tcond-TL)

        #aggregates. In this case gamma is taken from tensile stress
        #model from Gundlach et al. 2018
        elif self.flag_cond==2:
            a1=2.93E3
            a2=2.61E3
            alpha=2.73E-3
            
            sigma_cor_ref=(a1/(1-0.5)-a2)*0.5
            sigma_cor=(a1/(1.-self.Vfrac_ag)-a2)*self.Vfrac_ag
            xi=sigma_cor/sigma_cor_ref

            sigma_fit=alpha/self.r_mono_ag

            #Note that this needs to be 'converted' into a surface tension like factor if used elsewhere
            self.gamma=sigma_fit/xi#*self.Rcond/self.prefac

            

        else:
            print('ERROR \n Unknown flag_cond \n EXITING')
            sys.exit()


    ###################################################
    ###################################################
    #function to calculate the viscosity of the vapor
    def calc_visc(self):
        #H2, He: Parameters taken from fitting of experimental data
        if self.flag_vap==0:
            XH2=self.vap_prop[0]
            fit_H2=[ 0.01694574,  0.6938941 ]
            fit_He=[ 0.03818835,  0.69396091]

            #calculate the individual viscosity
            viscH2=fit_H2[0]*(self.Tvap**fit_H2[1])
            viscHe=fit_He[0]*(self.Tvap**fit_He[1])

            #comine viscosities linearly
            self.visc=(XH2*viscH2+(1-XH2)*viscHe)*1E-5

        #H2O: Taken from expression from Huber et al. 2009
        elif self.flag_vap==1:
            mustar=1E-6
            Tstar=647.096
            rhostar=322.0
            pstar=22.064E6

            Tbar=self.Tvap/Tstar
            rhobar=self.rho_vap/rhostar

            #mu0
            H0i=[1.67752,2.20462,0.6366564, -0.241605]
            temp=0
            for i in np.arange(4):
                temp+=H0i[i]/(Tbar**i)
            mu0=100.0*np.sqrt(Tbar)/temp

            #mu1
            H1ij=[[5.20094E-1,8.50895E-2,-1.08374,-2.89555E-1,0.0,0.0],
                  [2.22531E-1,9.99115E-1,1.88797,1.26613,0.0,1.20573E-1],
                  [-2.81378E-1,-9.06851E-1,-7.72479E-1,-4.89837E-1,-2.57040E-1,0.0],
                  [1.61913E-1,2.57399E-1,0.0,0.0,0.0,0.0],
                  [-3.25372E-2,0.0,0.0,6.98452E-2,0.0,0.0],
                  [0.0,0.0,0.0,0.0,8.72102E-3,0.0],
                  [0.0,0.0,0.0,-4.35673E-3,0.0,-5.93264E-4]]
            exponent=0.0
            for i in np.arange(6):
                temp=0.0
                for j in np.arange(7):
                    #print i, j, H1ij[j][i]
                    temp+=H1ij[j][i]*((rhobar-1)**j)
                exponent+=temp*(((1.0/Tbar)-1)**i)

            mu1=np.exp(rhobar*exponent)

            self.visc=mu0*mu1*mustar

        else:
            print('ERROR \n Unknown flag_vap \n EXITING')
            sys.exit()

    ###################################################
    #function to calculate the viscosity of the condensate
    def calc_visc_cond(self):


        if self.flag_cond==2: #if an aggregate assume inviscid
            self.visc_cond=0.0
            return
            

        #self.visc_cond=1E-2
        #return
    
        #Giordano et al. 2011 silicate EOS
        if self.flag_visc_cond_model==0:
            #parameters for the model
            A=-4.55
            bi=[159.6,-173.3,72.1,75.7,-39.,-84.1,141.5,-2.43,-0.91,17.6]
            ci=[2.75,15.7,8.3,10.2,-12.3,-99.5,0.3] 
            
            mmw={"SiO2": 28.09+2*16,
                 "TiO2": 47.88+2*16,
                 "Al2O3": 26.98*2+3*16,
                 "FeO": 55.85+16,
                 "MnO": 54.94+16,
                 "MgO": 24.31+16,
                 "CaO": 40.08+16,
                 "Na2O": 22.99*2+16,
                 "K2O": 39.1*2+16,
                 "P2O5": 30.97*2+5*16,
                 "H2O": 1.01*2+16,
                 "F2O-1": 19.*2-16}
            
            #Initialize an empty composition
            comp0={"SiO2": 0.0,
                   "TiO2": 0.0,
                   "Al2O3": 0.0,
                   "FeO": 0.0,
                   "MnO": 0.0,
                   "MgO": 0.0,
                   "CaO": 0.0,
                   "Na2O": 0.0,
                   "K2O": 0.0,
                   "P2O5": 0.0,
                   "H2O": 0.0,
                   "F2O-1": 0.0}
            
            #And merge in give comp
            comp0.update(self.comp)
            comp=comp0
    
            #renormalize the composition
            mtot=np.sum(list(comp.values()))
            for key in list(comp.keys()):
                comp[key]=comp[key]*100./mtot
                
            #If the composition given in wt% then convert to mol%
            if self.flag_comp_mol==0:
                for key in list(comp.keys()):
                    comp[key]=comp[key]/mmw[key]
                    
                #Renormalize to 100%
                Ntot=np.sum(list(comp.values()))
                for key in list(comp.keys()):
                    comp[key]=comp[key]*100./Ntot
                
        
            #Calculate the B term
            B=bi[0]*(comp["SiO2"]+comp["TiO2"])
            B+=bi[1]*comp["Al2O3"]
            B+=bi[2]*(comp["FeO"]+comp["MnO"]+comp["P2O5"])
            B+=bi[3]*comp["MgO"]
            B+=bi[4]*comp["CaO"]
            B+=bi[5]*(comp["Na2O"]+comp["H2O"]+comp["F2O-1"])
            B+=bi[6]*(comp["H2O"]+comp["F2O-1"]+np.log(1+comp["H2O"]))
            B+=bi[7]*(comp["SiO2"]+comp["TiO2"])*(comp["FeO"]+comp["MnO"]+comp["MgO"])
            B+=bi[8]*(comp["SiO2"]+comp["TiO2"]+comp["Al2O3"]+comp["P2O5"])*(comp["Na2O"]+comp["K2O"]+comp["H2O"])
            B+=bi[9]*comp["Al2O3"]*(comp["Na2O"]+comp["K2O"])
        
            #Calculate the C term
            C=ci[0]*comp["SiO2"]
            C+=ci[1]*(comp["TiO2"]+comp["Al2O3"])
            C+=ci[2]*(comp["FeO"]+comp["MnO"]+comp["MgO"])
            C+=ci[3]*comp["CaO"]
            C+=ci[4]*(comp["Na2O"]+comp["K2O"])
            C+=ci[5]*np.log(1+comp["H2O"]+comp["F2O-1"])
            C+=ci[6]*(comp["Al2O3"]+comp["FeO"]+comp["MnO"]+comp["MgO"]+comp["CaO"]-comp["P2O5"])\
                *(comp["Na2O"]+comp["K2O"]+comp["H2O"]+comp["F2O-1"])
            
            self.visc_cond=np.exp(A+B/(self.Tcond-C))
        
        else:
            raise ValueError("Condensate viscosity model unknown")
            self.visc_cond=np.nan

        return
    

    ###########################################################
    ###########################################################
    ###########################################################
    #functions for dynamical properties

    ###################################################
    ###################################################
    #function to calculate what regime the droplet is in
    def calc_drag_regime(self):
        self.MFP=const.k*self.Tvap/(np.sqrt(2)*self.sigma_coll*self.pvap)
        if (self.Rcond<((9.0/4.0)*self.MFP)):
            self.flag_drag_regime=0
        else:
            self.flag_drag_regime=1

    ###################################################
    ###################################################
    #function to calculate the Reynolds number
    def calc_Re(self,Rcond,vvap):
        Re=2.0*Rcond*self.rho_vap*(vvap)/self.visc
        return Re
    
    ###################################################
    ###################################################
    #function to calculate the drag coefficient. Does so depending on flag
    #0: Epstein regime. Equation 4 and 5 in Morris 2012
    #1: Low-v dynamic regime. Equation 19 in Brown & Lawler 2003
    #2: High-v dynamic regime. Equation 4 in Tanigawa et a. 2014
    #3: High-v dynamic regime. Equation 5 in Nagasawa et al. (on ArXiv). !!!! Probably incorrect!!!!
    #4: High-v dynamic regime. As used by Phil. Brasser+ 2007
    def calc_CD(self,Rcond,vvap,flag_CD_form=1):
        
        if flag_CD_form==0:
            s=vvap/np.sqrt(2.0*const.R*self.Tvap/(self.ma))
            CD=(2.0/(3.0*s))*np.sqrt(np.pi*self.Tcond/self.Tvap)\
                +(2.0*(s**2)+1.0)/(np.sqrt(np.pi)*(s**3))*np.exp(-(s**2))\
                +(4.0*(s**4)+4.0*(s**2)-1.0)/(2.0*(s**4))*math.erf(s)

            return CD

        elif (flag_CD_form==1): #include low velocity limit
            Re=self.calc_Re(Rcond,vvap)
            CD=24.0*(1.0+0.15*(Re**0.681))/Re + 0.407/(1.0+(8710.0/Re))
            return CD

        elif (flag_CD_form==2)|(flag_CD_form==3)|(flag_CD_form==4):

            #calculate the MAch number
            Mach=vvap/self.cs

            #calc Re and corresponding correction for flags 2 and 3
            Re=self.calc_Re(Rcond,vvap)
            correction=np.ones(np.size(Re))*0.2
            correction[np.where(Re<2E5)[0]]=0.4

            if flag_CD_form==2:
                CD=(((24.0/Re)+(40.0/(10.0+Re)))**(-1)+0.23*Mach)**(-1) + (2.0-correction)*Mach/(1.6+Mach)+correction
            elif flag_CD_form==3:
                CD=(((24.0/Re)+(40.0/(10.0+Re)))**(-1)+0.23*Mach)**(-2) + 2.0*(0.8*correction+Mach)/(1.6+Mach)
            elif flag_CD_form==4:

                #if the length of Re due to v_vap then go through each mach number
                if np.size(vvap)>1:
                    CD=zeros(np.size(Re))
                    temp=np.where(Mach>=1)[0]
                    if np.size(temp)>=1:
                        CD[temp]=2.0
                        temp=np.where((Mach<1)&(Re>=1E3))[0]
                    if np.size(temp)>=1:
                        CD[temp]=0.44+1.56*(Mach[temp]**2)
                        temp=np.where((Mach<1)&(Re<1E3))[0]
                    if np.size(temp)>=1:
                        CD[temp]=2*(Mach[temp]**2)+(1-Mach[temp]**2)*(24*(1+0.15*(Re[temp]**0.687)))/Re[temp]

                #if the length of Re is one
                elif np.size(Re)==1:
                    if Mach>=1:
                        CD=2.0
                    elif Re>=1E3:
                        CD=0.44+1.56*(Mach**2)
                    else:
                        CD=2*(Mach**2)+(1-Mach**2)*(24*(1+0.15*(Re**0.687)))/Re

                #if the length of Re due to any of the constituents of Re
                else:
                    if Mach>=1:
                        CD=2.0*np.ones(np.size(Re))
                    else:
                        CD=zeros(np.size(Re))
                        temp=np.where(Re>=1E3)[0]
                        if np.size(temp)>=1:
                            CD[temp]=0.44+1.56*(Mach**2)
                            temp=np.where(Re<1E3)[0]
                        if np.size(temp)>=1:
                            CD[temp]=2*(Mach**2)+(1-Mach**2)*(24*(1+0.15*(Re[temp]**0.687)))/Re[temp]


            return CD

        else:
            print("ERROR \n Unknown flag_CD_form \n EXITING")
            sys.exit()

    ###################################################
    ###################################################
    #function to calculate the critical drag coeff
    def calc_vcrit(self):
        if self.flag_cond==2: #need to convert tensile strength into a surface tension
            vcrit=2.0*self.prefac*(self.gamma*self.Rcond/self.prefac)/(24.*self.visc)
        else:
            vcrit=2.0*self.prefac*self.gamma/(24.*self.visc)
        self.vcrit=vcrit

    ####################################################
    ####################################################
    ##function to calculate the breakup radius based on Theofanous 2011
    #def calc_Rbreakup_Wec(self):
    #    self.Rbreakup_Wec=self.Wec*self.gamma/(self.rho_vap*(self.vvap**2))/2.0

    ####################################################
    ####################################################
    ##function to calculate the breakup radius based on Theofanous 2011
    #def calc_tbreakup_Wec(self):
    #    self.tbreakup=2.0*self.Rcond/self.vvap*np.sqrt(self.rho_cond/self.rho_vap) #using the prescription from Theo#fanous
    #    self.tbreakup_cs=2*self.Rcond/(self.cs_cond)

    def calc_Rbreakup_Wec_func(self,Rcond):

        We=self.rho_vap*self.vvap**2*2*Rcond/self.gamma
        Oh=self.visc_cond/np.sqrt(self.rho_cond*self.gamma*2*Rcond)

        return 12*(1+(4./3)*Oh**1.5)-We

    def calc_Rbreakup_Wec_inv_func(self,Rcond):

        We=self.rho_vap*self.vvap**2*2*Rcond/self.gamma
        Oh=0.0

        return 12*(1+(4./3)*Oh**1.5)-We

    def calc_Rbreakup_Wec2_func(self,Rcond):

        We=self.rho_vap*self.vvap**2*2*Rcond/self.gamma
        Oh=self.visc_cond/np.sqrt(self.rho_cond*self.gamma*2*Rcond)

        return 1E2+2.4E3*Oh-We

    def calc_breakup_Wec(self,vvap,Rcond):

        #calculate the critical We numbers
        self.Wec=12*(1+(4./3)*self.Oh**1.5)
        self.Wec2=1E2+2.4E3*self.Oh

        #Radius of particle at the critical We
        self.Rbreakup_Wec=self.Wec*self.gamma/(self.rho_vap*(vvap**2))/2.0

        sol=fsolve(self.calc_Rbreakup_Wec_func, self.Rbreakup_Wec, full_output=1,xtol=1E-14,maxfev=100000)
        self.Rbreakup_Wec=sol[0]

        sol=fsolve(self.calc_Rbreakup_Wec_inv_func, self.Rbreakup_Wec, full_output=1,xtol=1E-14,maxfev=100000)
        self.Rbreakup_Wec_inv=sol[0]

        sol=fsolve(self.calc_Rbreakup_Wec2_func, self.Rbreakup_Wec, full_output=1,xtol=1E-14,maxfev=100000)
        self.Rbreakup_Wec2=sol[0]

        #the sound speed breakup time
        self.tbreakup_cs=2*Rcond/(self.cs_cond)

        #determine if the droplet is either in the Epstein regime or stable
        if self.flag_drag_regime==0:
            self.flag_breakup_regime=-1
            self.tbreakup=self.tbreakup_cs
            return

        elif (self.flag_drag_regime==1)&(self.We<=self.Wec): #if stable this is really easy
            self.flag_breakup_regime=0
            #self.tbreakup_cs=np.inf
            self.tbreakup=np.inf

            return


        

        #find the RTP breakup time
        if (self.Oh<1.):
            tbreakup_RTP=self.tI_plus
            self.flag_breakup_regime=1
        else:
            tbreakup_RTP=self.regimeII_prefac*self.tv_plus
            self.flag_breakup_regime=4

        #now find what domain we are in
        if (self.We>self.Wec)&(self.We<self.Wec2): #Regime I or IV

            self.tbreakup=tbreakup_RTP
            #leave the flag_breakup_regime the same

        elif (self.We>self.Wec2)&(self.Oh<0.1)&(self.We<1.5E3): #Regime II

            self.flag_breakup_regime=2

            self.tbreakup=np.amin([tbreakup_RTP,2*self.tI_plus])

        elif (self.We>self.Wec2)&(self.Oh<0.1): #Regime III

            self.flag_breakup_regime=3

            self.tbreakup=np.amin([tbreakup_RTP,0.6*self.tI_plus])

        elif (self.We>self.Wec2)&(self.Oh>0.1): #Regimes V and VI

            self.flag_breakup_regime=5

            self.tbreakup=np.amin([tbreakup_RTP,2.0*self.tI_plus])

        else:
            print(self.Rcond, self.Rbreakup_Wec, self.We, self.Wec, self.Wec2, self.Oh>0.1)
            raise ValueError("Unknown breakup time regime")

        return

        

            
        

    ###################################################
    ###################################################
    #function to calculate to solve for the equilibrium droplet size

    ############################################
    #we need this cost function
    def calc_Rcond_func(self,Rcond,vvap,flag_CD_form=1):
        CD=self.calc_CD(Rcond,vvap,flag_CD_form=flag_CD_form)
        Rcond_calc=(self.prefac*self.gamma/self.rho_vap/(vvap**2)/CD)
        return (Rcond_calc-Rcond)/1E-3

    def calc_Rcond_func_ag(self,Rcond,vvap,flag_CD_form=1):
        CD=self.calc_CD(Rcond,vvap,flag_CD_form=flag_CD_form)

        Fdrag=0.5*self.rho_vap*(vvap**2)*CD*np.pi*(Rcond**2)
        #note that gamma needs to be converted back to surface tension
        # Ftense=(self.gamma/self.Rcond*self.prefac)*np.pi*(Rcond**2)
        Ftense=(self.gamma)*np.pi*(Rcond**2)
        
        return (Fdrag-Ftense)#/Ftense

    ############################################
    #function to calculate the equilibrium droplet size
    def calc_Req(self,vvap,flag_CD_form=1, flag_method=0, xtoll=1E-14):

        #aggregates
        if self.flag_cond==2:
            if flag_CD_form==0:
                CD=self.calc_CD(0.0,vvap,flag_CD_form=0)

                Fdrag=0.5*self.rho_vap*(vvap**2)*CD*np.pi*(self.Rcond**2)
                #note that gamma could need to be converted back to surface tension
                # Ftense=(self.gamma/self.Rcond*self.prefac)*np.pi*(self.Rcond**2)
                Ftense=(self.gamma)*np.pi*(self.Rcond**2)

                
                if Fdrag<Ftense:
                    self.flag_FB_stable=True
                    Req=self.Rcond
                else:
                    self.flag_FB_stable=False
                    Req=np.nan
                    
                
                return Req, 0

            else:
                # print('Aggregates not in the epstein')
                #else find eq droplet size
                count=0
                Ncount=10
                Rinit=np.logspace(-8,2,Ncount)
                while count<Ncount:
                    sol=fsolve(self.calc_Rcond_func_ag, Rinit[count], args=(vvap,flag_CD_form),full_output=1,xtol=xtoll,maxfev=100000)
                    Req=sol[0]
                    #print(Req,self.calc_Rcond_func(Req,vvap,flag_CD_form))
                    if (np.isnan(Req))|(Req<=0.0)|(sol[2]!=1)|(self.calc_Rcond_func_ag(Req,vvap,flag_CD_form=flag_CD_form)>xtoll*100):
                        Req=np.nan
                        count+=1
                    else:
                        break

                if count==Ncount:
                    # print('Satisfactory solution not found',self.calc_Rcond_func_ag(Req,vvap,flag_CD_form=flag_CD_form))
                    return np.nan, -1

                if self.Rcond<Req:
                    self.flag_FB_stable=True
                else:
                    self.flag_FB_stable=False

                return Req, 0

            
        #spheres/droplets
        #if in Epstein regime then the solution is simple
        if flag_CD_form==0:
            CD=self.calc_CD(0.0,vvap,flag_CD_form=0)
            Req=self.prefac*self.gamma/(CD*self.rho_vap*(vvap**2))
            if self.Rcond<Req:
                self.flag_FB_stable=True
            else:
                self.flag_FB_stable=False
            return Req, 0

        #if in sub-mac dynamic regime test to see if above critical velocity. If so, return nan
        if (flag_CD_form==1)&(vvap>=self.vcrit):
            Req=np.nan
            flag_crit=1
            self.flag_FB_stable=False
        else:
            #else find eq droplet size
            flag_crit=0
            count=0
            Ncount=100
            Rinit=np.logspace(-4,4,Ncount)
            while count<Ncount:
                sol=fsolve(self.calc_Rcond_func, Rinit[count], args=(vvap,flag_CD_form),full_output=1,xtol=xtoll,maxfev=10000)
                Req=sol[0]
                #print(Req,self.calc_Rcond_func(Req,vvap,flag_CD_form))
                if (np.isnan(Req))|(Req<=0.0)|(sol[2]!=1)|(self.calc_Rcond_func(Req,vvap,flag_CD_form=flag_CD_form)>(xtoll*100)):
                    Req=np.nan
                    count+=1
                else:
                    break
                
            if self.Rcond<Req:
                self.flag_FB_stable=True
            else:
                self.flag_FB_stable=False
            
        return Req,flag_crit

    ###################################################
    ###################################################
    #function to solve for when Epstein and dynamic drag are the same

    ############################################
    #we need this cost function
    def calc_regime_crossover_func(self,vvap,flag_CD_form=1):
        Requib_Ep=self.calc_Req(vvap,flag_CD_form=0)[0]
        temp=self.calc_Req(vvap,flag_CD_form=flag_CD_form)
        if temp[1]==1:
            return (self.vcrit-vvap)-10
        else:
            Requib_dy=temp[0]
            return (np.log10(Requib_dy)-np.log10(Requib_Ep))

    ############################################
    def calc_regime_crossover(self,flag_CD_form=1):
        sol=newton(self.calc_regime_crossover_func, 0.1E3, args=(flag_CD_form,),maxiter=1000)
        vcross=sol
        Requib_cross=self.calc_Req(vcross,flag_CD_form=0)[0]
        return vcross, Requib_cross

    
    ###########################################################
    ###########################################################
    ###########################################################
    #functions to calculate a time evolving body

    ###################################################
    ###################################################
    #function that calculates the velocity evolution assuming the radius stays constant

    ############################################
    #function that calculates the time derivative assuming the radius stays constant
    #z=[drdt, r]
    def calc_time_derivative_constR(self, z, t, vvap_func,vvap_func_params,Rcond,flag_CD_form=1,flag_stab_check=0):
        #unpack z for easy of use
        vdrop=z[0]
        xdrop=z[1]
        vvap=z[2]
        xvap=z[3]

        #initialize array
        dzdt=np.empty(4)

        #acceleration of droplet given by
        if vvap==vdrop:
            CD=0.0
            dzdt[0]=0.0
        else:
            
#             test=self.flag_drag_regime
#             self.Rcond=Rcond
#             # self.vvap=vvap-vdrop
#             self.calc_drag_regime()
            
#             if self.flag_drag_regime!=test:
#                 print('We have swapped drage regime')
        
#             if self.flag_CD_self_consistent==1:
#                 if self.flag_drag_regime==0:
#                     flag_CD_form=0
#                     self.flag_CD_form=0
#                 else:
#                     flag_CD_form=self.flag_CD_form_dynamic
#                     self.flag_CD_form=self.flag_CD_form_dynamic
#             else:
#                 flag_CD_form=self.flag_CD_form_dynamic
#                 self.flag_CD_form=self.flag_CD_form_dynamic
                
                
            CD=self.calc_CD(Rcond,np.absolute(vvap-vdrop),flag_CD_form=flag_CD_form)
            #dzdt[0]=np.sign(vvap-vdrop)*0.5*CD*np.pi*self.rho_vap*((Rcond*(vvap-vdrop))**2)/self.m

            dzdt[0]=np.sign(vvap-vdrop)*1.5*CD*self.rho_vap*(((vvap-vdrop))**2)/(4.0*self.rho_cond*(Rcond))
            
            if flag_stab_check>0:
                if self.flag_stop_stab==1:
                    if (self.flag_drag_regime==0)|(self.flag_cond==2)|(flag_stab_check==1):
                        temp=self.calc_Req(vvap,flag_CD_form=self.flag_CD_form)
                        if (temp[0]<Rcond)|(temp[1]!=0):
                            self.flag_stop_stab=0

                    else:
                        self.We=self.rho_vap*(vvap-vdrop)**2*2*Rcond/self.gamma
                        self.Re=self.calc_Re(Rcond,vvap-vdrop)
                        self.Oh=self.visc_cond/np.sqrt(self.rho_cond*self.gamma*2*Rcond)
                        self.calc_breakup_Wec(vvap,Rcond)

                        if (self.flag_breakup_regime!=0):
                            self.flag_stop_stab=0

        #the time derivative of x is the same as 
        dzdt[1]=vdrop

        #find the vapor velocity
        temp=vvap_func(t,z,vvap_func_params)
        dzdt[2]=temp[0]
        dzdt[3]=vvap

        #print('\t\t\t',CD,z)

        return dzdt

    ############################################
    #function for acceleratioin of vapor
    #assume constant velocity gas
    def const_vvel(self,t,z,vvap_func_params):
        return 0.0, 0

    ############################################
    #function for acceleratioin of vapor
    #assume constant deceleration
    def const_vacc(self,t,z,vvap_func_params):
        return vvap_func_params[0], 0


    ############################################
    #function that calculates the time derivative assuming the radius stays constant
    def calc_path_constR(self,tcalc,zinit,vvap_func,vvap_func_params,Rcond,flag_CD_form=1,flag_stab_check=0):
        
        #solve the odeint
        sol = odeint(self.calc_time_derivative_constR, zinit, tcalc, args=(vvap_func,vvap_func_params,Rcond,flag_CD_form,flag_stab_check),rtol=1E-10)

            
        return sol

    ###################################################
    ###################################################
    #function that calculates the point at which a particle reaches a certain velocity
    #evolution assuming the radius stays constant
    def find_stopping_time(self,vstop,tguess,toll,Ncalc,Ncalc_final,flag_CD_form=1,flag_stab_check=0):
        #initialize the test and stop times
        test=1.0
        tbegin=1E-15
        tstop=tguess
        count=0
        
        while test>toll:
            #integrate the function
            tcalc=np.append([0.0],np.linspace(tbegin,tstop,Ncalc))
            #print(tcalc)
            #print(tguess)
            #print(tcalc)
            
            sol=self.calc_path_constR(tcalc,[self.vvap,0.0,0.0,0.0],self.const_vvel,\
                                      [],self.Rcond,flag_CD_form=flag_CD_form,flag_stab_check=flag_stab_check)
            #if the function didn't reach the desired velocity then increase tstop
            if sol[-1][0]>vstop:
                Ncalc+=1 #stops the function getting in a loop
                #print('that')
                tbegin=np.maximum(1E-15,tbegin*(1-(toll/2.0/vstop)))#scale to vstop to get dimensions sort of right
                #print('\t\t\t',tstop, np.amin([1.0,np.absolute(test)]), test)
                tstop=tstop*(1+np.amin([1.0,np.absolute(test)]))
                test=np.absolute((sol[-1][0]-vstop)/vstop)
                #print(tbegin,tstop,tbegin-tstop)
            else:
                if test<toll*100: #increase the step size to increase resolution near the stopping point
                    #Ncalc=np.maximum(200,int(1.0/toll)*10)
                    Ncalc+=100*int(np.absolute(np.log10(toll)))+1
                else:
                    #Ncalc=np.maximum(200,int(1.0/test))
                    Ncalc+=100*int(np.absolute(np.log10(test)))+1
                
                #Ncalc+=501 #chosen roughly to give the fastest convergence
                temp=np.where(sol[:,0]>vstop)[0]
                if np.size(temp)<2:
                    tbegin=np.maximum(1E-15,tbegin*(1-(toll/2.0/vstop))) #scale to vstop to get dimensions sort of right
                    tstop=tcalc[temp[-1]]
                else:
                    tbegin=tcalc[temp[-2]]
                    tstop=tcalc[temp[-1]]
                #print(np.size(temp),temp[-1],sol[0:10,0])
                #print(tbegin,tstop,tbegin-tstop)
                test=np.absolute((sol[temp[-1],0]-vstop)/vstop)
            count+=1
        #print(count,tstop,test)
                
            #print(test,Ncalc,tbegin,tstop)
        #print('\t',tstop)
        #record properties
        self.toll_stop=toll
        self.Ncalc_stop=Ncalc_final
        self.tstop=tstop
        self.tcalc_stop=np.linspace(0,tstop,self.Ncalc_stop)
        self.sol_stop=self.calc_path_constR(self.tcalc_stop,[self.vvap,0.0,0.0,0.0],self.const_vvel,\
                                            [],self.Rcond,flag_CD_form=flag_CD_form,flag_stab_check=flag_stab_check)
        

        #print(count)
        return(tstop)#print(tstop)


    ######################################################################
    # a new set of funcations to find stopping time that are more computationally efficient/robust in most cases
    ############################################
    #function that calculates the time derivative assuming the radius stays constant
    #z=[drdt, r]
    def calc_time_derivative_constR2(self, t, z, vvap_func,vvap_func_params,Rcond,vstop,flag_CD_form=1,flag_stab_check=0):
        
        
        #unpack z for easy of use
        vdrop=z[0]
        xdrop=z[1]
        vvap=z[2]
        xvap=z[3]

        #initialize array
        dzdt=np.empty(4)

        #acceleration of droplet given by
        if vvap==vdrop:
            CD=0.0
            dzdt[0]=0.0
        else:
#             test=self.flag_drag_regime
#             self.Rcond=Rcond
#             self.calc_drag_regime()
            
#             if self.flag_drag_regime!=test:
#                 print('We have swapped drage regime')
        
#             if self.flag_CD_self_consistent==1:
#                 if self.flag_drag_regime==0:
#                     flag_CD_form=0
#                     self.flag_CD_form=0
#                 else:
#                     flag_CD_form=self.flag_CD_form_dynamic
#                     self.flag_CD_form=self.flag_CD_form_dynamic
#             else:
#                 flag_CD_form=self.flag_CD_form_dynamic
#                 self.flag_CD_form=self.flag_CD_form_dynamic
            
            CD=self.calc_CD(Rcond,np.absolute(vvap-vdrop),flag_CD_form=flag_CD_form)
            #dzdt[0]=np.sign(vvap-vdrop)*0.5*CD*np.pi*self.rho_vap*((Rcond*(vvap-vdrop))**2)/self.m

            dzdt[0]=np.sign(vvap-vdrop)*1.5*CD*self.rho_vap*(((vvap-vdrop))**2)/(4.0*self.rho_cond*(Rcond))
            
            
           
        #the time derivative of x is the same as 
        dzdt[1]=vdrop

        #find the vapor velocity
        temp=vvap_func(t,z,vvap_func_params)
        dzdt[2]=temp[0]
        dzdt[3]=vvap

        #print('\t\t\t',CD,z)

        return dzdt

   
    ###############################################
    def find_stopping_time2(self,vstop,tguess,toll,Ncalc,Ncalc_final,flag_CD_form=1,vvap_func=None,vvap_func_params=[],flag_stab_check=0):

        #need this function for finding when the particle stops
        def event_stopped(t,z, vvap_func,vvap_func_params,Rcond,vstop,flag_CD_form=1,flag_stab_check=flag_stab_check):
            return z[0]-vstop
        event_stopped.terminal = True
        
        if vvap_func==None:
            vvap_func=self.const_vvel

        #Need to have a range of t that includes the stopping time
        success=False
        count=0
        while success==False:
            tspan=[0.,tguess]
        
            sol=solve_ivp(self.calc_time_derivative_constR2, tspan, [self.vvap,0.0,0.0,0.0], args=(vvap_func,vvap_func_params,self.Rcond,vstop,self.flag_CD_form, flag_stab_check),rtol=1E-10,atol=1E-10,events=event_stopped,method='LSODA')#,t_eval=np.linspace(0,1200,1000000),dense_output=True)

            #If not hit the event then increase the max t, scaling based on a roughly log-linear v-t relationship
            if np.size(sol.t_events)==0:
                tguess=np.amax([1.5*tguess,tguess+sol.y[1,-1]/(np.log10(self.vvap)-np.log10(sol.y[0,-1]))])
                #print(self.Rcond,tguess)
            else:
                success=True


        #record the results
        self.tstop=sol.t_events[0][0]
        self.Lstop=sol.y_events[0][0][1]
        
        self.sol_stop=np.transpose(sol.y)
        self.tcalc_stop=sol.t
        
        #test stability
        self.test_stop_stab(flag_stab_check)

        return(sol.t_events[0][0])
    
    
    ###############################################
    def find_stopping_time3(self,vstop,tguess,toll,Ncalc,Ncalc_final,flag_CD_form=1,vvap_func=None,vvap_func_params=[],flag_stab_check=0):

        #need this function for finding when the particle stops
        #need this function for finding when the particle stops
        def event_stopped(t,z, vvap_func,vvap_func_params,Rcond,vstop,flag_CD_form=1,flag_stab_check=flag_stab_check):
            return z[2]-z[0]-vstop
        event_stopped.terminal = True
        
        if vvap_func==None:
            vvap_func=self.const_vvel

        #Need to have a range of t that includes the stopping time
        success=False
        count=0
        while success==False:
            tspan=[0.,tguess]
        
            sol=solve_ivp(self.calc_time_derivative_constR2, tspan, [0.0,0.0,self.vvap,0.0], args=(vvap_func,vvap_func_params,self.Rcond,vstop,self.flag_CD_form, flag_stab_check),rtol=1E-10,atol=1E-10,events=event_stopped,method='DOP853')#,t_eval=np.linspace(0,1200,1000000),dense_output=True)

            #If not hit the event then increase the max t, scaling based on a roughly log-linear v-t relationship
            if np.size(sol.t_events)==0:
                tguess=np.amax([1.5*tguess,tguess+sol.y[1,-1]/(np.log10(self.vvap)-np.log10(sol.y[0,-1]))])
                #print(self.Rcond,tguess)
            else:
                success=True


        #record the results
        self.tstop=sol.t_events[0][0]
        self.Lstop=sol.y_events[0][0][3]-sol.y_events[0][0][1]
        
        self.sol_stop=np.transpose(sol.y)
        self.tcalc_stop=sol.t
        
        #test stability
        self.test_stop_stab(flag_stab_check)

        return(sol.t_events[0][0])
    
    
    
    
    ###############################################
    def find_stopping_time_2step(self,vstop,tguess,toll,Ncalc,Ncalc_final,flag_CD_form=1,vvap_func=None,vvap_func_params=[],deltav_L=None,flag_vap2=None,rho_vap2=None,Tvap2=None,deltav2=None, flag_stab_check=0):

        #if not passing second wave information, just take the properties of the first wave or leave as set values
        if (flag_vap2==None)&(self.flag_vap2==None):
            self.flag_vap2=self.flag_vap
        elif (flag_vap2!=None):
            self.flag_vap2=flag_vap2
            
        if (rho_vap2==None)&(self.rho_vap2==None):
            self.rho_vap2=self.rho_vap
        elif (rho_vap2!=None):
            self.rho_vap2=rho_vap2
            
        if (Tvap2==None)&(self.Tvap2==None):
            self.Tvap2=self.Tvap
        elif (Tvap2!=None):
            self.Tvap2=Tvap2
            
        if (deltav_L==None)&(self.deltav_L==None):
            self.deltav_L=0.0
        elif (deltav_L!=None):
            self.deltav_L=deltav_L
            
        if (deltav2==None)&(self.deltav2==None):
            self.deltav2=0.0
        elif (deltav2!=None):
            self.deltav2=deltav2
            
        vvap_func_params=np.append(vvap_func_params,[self.deltav_L])
        
        #need this function for finding when the particle reaches the second wave
        def event_step(t,z, vvap_func,vvap_func_params,Rcond,vstop,flag_CD_form=1,flag_stab_check=flag_stab_check):
            return z[1]-z[3]+vvap_func_params[-1]
        event_step.terminal = True
        
        #need this function for finding when the particle stops
        def event_stopped(t,z, vvap_func,vvap_func_params,Rcond,vstop,flag_CD_form=1,flag_stab_check=flag_stab_check):
            return z[2]-z[0]-vstop
        event_stopped.terminal = True
        
        if vvap_func==None:
            vvap_func=self.const_vvel

        #Need to have a range of t that includes the stopping time
        success=False
        count=0
        while success==False:
            tspan=[0.,tguess]
        
            sol=solve_ivp(self.calc_time_derivative_constR2, tspan, [0.0,0.0,self.vvap,0.0], args=(vvap_func,vvap_func_params,self.Rcond,vstop,self.flag_CD_form, flag_stab_check),rtol=1E-10,atol=1E-10,events=(event_stopped,event_step),method='DOP853')#,t_eval=np.linspace(0,1200,1000000),dense_output=True)
            
            # print(sol)

            #If not hit the event then increase the max t, scaling based on a roughly log-linear v-t relationship
            if (np.size(sol.t_events[0])==0)&(np.size(sol.t_events[1])==0):
                tguess=np.amax([1.5*tguess,tguess+sol.y[1,-1]/(np.log10(self.vvap)-np.log10(sol.y[0,-1]))])
                #print(self.Rcond,tguess)
            else:
                success=True
                
            

        # print(sol.y)
        if (np.size(sol.t_events[1])==0): #if not reached second shock
            # print('coupled before kick')
            self.flag_hit_2nd_wave=0
            
            self.tstop1=sol.t_events[0][0]
            # self.Lstop1=sol.y_events[0][0][1]
            #put in ref frame relative to first wave
            self.Lstop1=sol.y_events[0][0][3]-sol.y_events[0][0][1]
            
            self.tstop2=0.0
            self.Lstop2=0.0
           
            self.tstop=self.tstop1
            self.Lstop=self.Lstop1

            self.sol_stop=np.transpose(sol.y)
            self.tcalc_stop=sol.t

        else:
            
            self.flag_hit_2nd_wave=1
            
            self.tstop1=sol.t_events[1][0]
            # self.Lstop1=sol.y_events[1][0][1]
            #put in ref frame relative to first wave
            self.Lstop1=sol.y_events[1][0][3]-sol.y_events[1][0][1]

            self.sol_stop=np.transpose(sol.y)
            self.tcalc_stop=sol.t
            
            #record the initial properties of the gas
            (vvap_save,flag_vap_save,rho_vap_save,Tvap_save,Tcond_save)=(self.vvap, \
                                                                         self.flag_vap,self.rho_vap,self.Tvap,self.Tcond)
            
            #reintialize droplet
            self.init_droplet(self.flag_cond, self.Rcond,self.Tvap2,self.flag_vap2,self.vap_prop, self.rho_vap2,self.Tvap2,self.deltav2,self.cs, \
                            self.prefac,self.flag_CD_form_dynamic,\
                              deltav_L=self.deltav_L, deltav2=self.deltav2, flag_vap2=self.flag_vap2, \
                              rho_vap2=self.rho_vap2, Tvap2=self.Tvap2,\
                             Vfrac_ag=self.Vfrac_ag, rho_mono_ag=self.rho_mono_ag, r_mono_ag=self.r_mono_ag)
            
            
            #Need to have a range of t that includes the stopping time
            success=False
            count=0
            start_t=sol.t_events[1][0]
            start_y=sol.y_events[1][0]
            start_y[2]+=self.deltav2
            tguess+=sol.t_events[1][0]
            while success==False:
                tspan=[start_t,tguess]

                sol=solve_ivp(self.calc_time_derivative_constR2, tspan, start_y, args=(vvap_func,vvap_func_params,self.Rcond,vstop,self.flag_CD_form, flag_stab_check),rtol=1E-10,atol=1E-10,events=(event_stopped),method='DOP853')#,t_eval=np.linspace(0,1200,1000000),dense_output=True)


                #If not hit the event then increase the max t, scaling based on a roughly log-linear v-t relationship
                if (np.size(sol.t_events[0])==0):
                    tguess=np.amax([1.5*tguess,tguess+sol.y[1,-1]/(np.log10(self.vvap)-np.log10(sol.y[0,-1]))])
                    #print(self.Rcond,tguess)
                else:
                    success=True
            

            #record the results
            self.tstop=sol.t_events[0][0]
            #put in ref frame relative to first wave
            self.Lstop=vvap_save*self.tstop-sol.y_events[0][0][1]
            
            #second start relative to first
            self.tstop2=sol.t_events[0][0]-self.tstop1
            #stop difference relative to second wave
            # self.Lstop2=sol.y_events[0][0][2]*(self.tstop-self.tstop1)-(sol.y_events[0][0][1]-self.sol_stop[-1,1])
            self.Lstop2=sol.y_events[0][0][2]*(self.tstop-self.tstop1)-(sol.y_events[0][0][1]-self.sol_stop[-1,1])


            # print(self.tcalc_stop)
            self.sol_stop=np.concatenate((self.sol_stop,np.transpose(sol.y)))
            self.tcalc_stop=np.append( self.tcalc_stop,sol.t)
            
            #return the original values for the vapor properties
            self.init_droplet(self.flag_cond, self.Rcond,Tcond_save,flag_vap_save,self.vap_prop, rho_vap_save,Tvap_save,vvap_save,self.cs, \
                              self.prefac,self.flag_CD_form_dynamic,\
                             deltav_L=self.deltav_L,deltav2=self.deltav2,flag_vap2=self.flag_vap2,rho_vap2=self.rho_vap2,Tvap2=self.Tvap2,\
                             Vfrac_ag=self.Vfrac_ag, rho_mono_ag=self.rho_mono_ag, r_mono_ag=self.r_mono_ag)
        
        #test the stability
        self.test_stop_stab(flag_stab_check)

        return(sol.t_events[0][0])

    
    #function to test stability of droplets as they strop
    def test_stop_stab(self,flag_stab_check):
        
        #record the initial properties of the gas
        (vvap_save,flag_vap_save,rho_vap_save,Tvap_save,Tcond_save)=(self.vvap, \
                                                                        self.flag_vap,self.rho_vap,self.Tvap,self.Tcond)
        
        self.flag_stop_stab=1
        for i in np.arange(np.size(self.tcalc_stop)):
            deltav=np.abs(self.sol_stop[i,0]-self.sol_stop[i,2])
            
            #if we are into the plume
            if (self.flag_hit_2nd_wave==1)&(self.tcalc_stop[i]>self.tstop1):
                #reintialize droplet
                self.init_droplet(self.flag_cond, self.Rcond,self.Tvap2,self.flag_vap2,self.vap_prop, self.rho_vap2,self.Tvap2,self.deltav2,self.cs, \
                                self.prefac,self.flag_CD_form_dynamic,\
                                  deltav_L=self.deltav_L, deltav2=self.deltav2, flag_vap2=self.flag_vap2, \
                                  rho_vap2=self.rho_vap2, Tvap2=self.Tvap2,\
                                 Vfrac_ag=self.Vfrac_ag, rho_mono_ag=self.rho_mono_ag, r_mono_ag=self.r_mono_ag)
            else:
                self.init_droplet(self.flag_cond, self.Rcond,Tcond_save,flag_vap_save,self.vap_prop, rho_vap_save,Tvap_save,vvap_save,self.cs, \
                    self.prefac,self.flag_CD_form_dynamic,\
                    deltav_L=self.deltav_L,deltav2=self.deltav2,flag_vap2=self.flag_vap2,rho_vap2=self.rho_vap2,Tvap2=self.Tvap2,\
                                 Vfrac_ag=self.Vfrac_ag, rho_mono_ag=self.rho_mono_ag, r_mono_ag=self.r_mono_ag)
                

            if (self.flag_drag_regime==0)|(self.flag_cond==2)|(flag_stab_check==1):
                temp=self.calc_Req(deltav,flag_CD_form=self.flag_CD_form)
                if (temp[0]<self.Rcond)|(temp[1]!=0):
                    self.flag_stop_stab=0
                    break


            else:
                self.We=self.rho_vap*(deltav)**2*2*self.Rcond/self.gamma
                self.Re=self.calc_Re(self.Rcond,deltav)
                self.Oh=self.visc_cond/np.sqrt(self.rho_cond*self.gamma*2*self.Rcond)
                self.calc_breakup_Wec(deltav,self.Rcond)

                if (self.flag_breakup_regime!=0):
                    self.flag_stop_stab=0
                    break
                    
                    
        #reestablish the droplet as was
        self.init_droplet(self.flag_cond, self.Rcond,Tcond_save,flag_vap_save,self.vap_prop, rho_vap_save,Tvap_save,vvap_save,self.cs, \
                  self.prefac,self.flag_CD_form_dynamic,\
                 deltav_L=self.deltav_L,deltav2=self.deltav2,flag_vap2=self.flag_vap2,rho_vap2=self.rho_vap2,Tvap2=self.Tvap2,\
                         Vfrac_ag=self.Vfrac_ag, rho_mono_ag=self.rho_mono_ag, r_mono_ag=self.r_mono_ag)
        
            
        return
                    