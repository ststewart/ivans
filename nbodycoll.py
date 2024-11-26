import numpy as np

def readcolllog(fname='ss.coll.txt'):
    """
    readcolllog(<fname=file>)
    - Read standard ss collision log
    - returns a tuple of arrays
    """
    t=[]
    id1=[]
    id2=[]
    m1=[]
    m2=[]
    x1=[]
    x2=[]
    y1=[]
    y2=[]
    z1=[]
    z2=[]
    vx1=[]
    vx2=[]
    vy1=[]
    vy2=[]
    vz1=[]
    vz2=[]
    mlr = []
    mslr = []
    nfrag = []
    mfrag = []
    m2dust = []
    with open(fname) as fo:
        lines = fo.readlines()
        for j in range(len(lines)):
            if lines[j].count('PLANETESIMAL')>1:
                mf = 0
                nf = 0
                mf0 = 0
                mf1 = 0
                mdeb = 0
                t.append(float(lines[j].split('=')[1]))
                d1 = lines[j+1].replace('=',',').split(',')
                d2 = lines[j+2].replace('=',',').split(',')
                id1.append(float(d1[7]))
                m1.append(float(d1[9]))
                x1.append(float(d1[17].strip('(')))
                y1.append(float(d1[18]))
                z1.append(float(d1[19].strip(')')))
                vx1.append(float(d1[21].strip('(')))
                vy1.append(float(d1[22]))
                vz1.append(float(d1[23].strip(')')))
                id2.append(float(d2[7]))
                m2.append(float(d2[9]))
                x2.append(float(d2[17].strip('(')))
                y2.append(float(d2[18]))
                z2.append(float(d2[19].strip(')')))
                vx2.append(float(d2[21].strip('(')))
                vy2.append(float(d2[22]))
                vz2.append(float(d2[23].strip(')')))
                d3 = lines[j+3]
                if d3.count('MERGE') + d3.count('BOUNCE') + d3.count('FRAG') >0:
                    d4 = lines[j+4].replace('=',',').split(',')
                    mf0 = float(d4[9])
                    mf += mf0
                    nf += 1
                if d3.count('BOUNCE') + d3.count('FRAG') >0:
                    d5 = lines[j+5].replace('=',',').split(',')
                    mf1 = float(d5[9])
                    mf += mf1
                    nf += 1
                if d3.count('FRAG') >0:
                    dd = lines[j+3+nf].replace('=',',').split(',')
                    while dd.count('***out')>0:
                        mf += float(dd[9])
                        nf += 1
                mlr.append(mf0)
                mslr.append(mf1)
                mfrag.append(mf)
                nfrag.append(nf)
                if abs((float(d1[9]) + float(d2[9]) - mf)/(float(d1[9]) + float(d2[9])))>1e-5:
                    mdeb = float(d1[9]) + float(d2[9]) - mf
                m2dust.append(mdeb)
            else:
                continue
    return np.array(t),np.array(id1),np.array(id2),np.array(m1),np.array(m2),np.array(x1),np.array(x2),np.array(y1),np.array(y2),np.array(z1),np.array(z2),np.array(vx1),np.array(vx2),np.array(vy1),np.array(vy2),np.array(vz1),np.array(vz2),np.array(mlr),np.array(mslr),np.array(nfrag),np.array(mfrag),np.array(m2dust)