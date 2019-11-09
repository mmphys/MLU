import numpy as np
import sys

#quarks=['l','h1','h2']
quarkL=['h2']
quarkR=['l']
quark3=['l'] #This one should always be light (in our case)

OldStyle=False
conf=3000
if OldStyle:
    # Older runs, where Felix made the unsmeared sinks, i.e. 3040, 3080, 3120
    PrefixUnsmeared  = '/home/dp008/dp008/dc-erbe1/run/heavy_light/output/MesonSink_read_'
    DirectorySmeared = '/home/dp008/dp008/dc-mars3/data/201907HL/smeared/'
else:
    # Newer runs, where I made everything, i.e. 3000, 3160
    PrefixUnsmeared  = '/home/dp008/dp008/dc-mars3/data/201910Plat/unsmeared/C1/MesonSink_read_'
    DirectorySmeared = '/home/dp008/dp008/dc-mars3/data/201910Plat/mesons/C1/Distil/'
print 'PrefixUnsmeared:  ' + PrefixUnsmeared
print 'DirectorySmeared: ' + DirectorySmeared
print

#Which momentum?
for p2 in range(0,5):
    nMom=0
    moms=[]
    moms_neg=[]
    for x in range(-p2,p2+1):
        for y in range(-p2,p2+1):
            for z in range(-p2,p2+1):
                if(x**2+y**2+z**2==p2):
                    nMom+=1
                    moms.append(str(x)+'_'+str(y)+'_'+str(z))
                    moms_neg.append(str(-x)+'_'+str(-y)+'_'+str(-z))

    print 'Momentum: ' + str(p2)
    print 'Count:    ' + str(nMom)
    print 'Positive: ' + str(moms)
    print 'Negative: ' + str(moms_neg)
    print

    #sys.exit(1)
    #continue

    # Which type of contraction?
    for ContractFor in range(0,3):
        if p2 != 2 and ContractFor != 0:
            continue
        if ContractFor >= 2:
            ContractFor = 2
            ContractPrefix = 'c2p0'
        elif ContractFor == 1:
            ContractPrefix = 'c3p0'
        else:
            ContractFor = 0
            ContractPrefix = 'c30p'

        BaseDir = ContractPrefix
        BaseDir += '_'
        BaseDir += quarkL[0] #Oh dear. I'm mixing lists and single strings ... bad
        BaseDir += '_'
        BaseDir += quarkR[0]
        BaseDir += '_p2_'
        BaseDir += str(p2)

        FileName  = BaseDir + '.xml'
        OutputDir = BaseDir
        TmpDir    = 'tmp_' + BaseDir

        print 'Creating XML: ' + FileName
        print 'OutputDir:    ' + OutputDir
        print 'TmpDir:       ' + TmpDir
        print

        indent=''
        out=open(FileName,'w')
        out.write('<?xml version="1.0"?>\n')
        out.write(indent + '<grid>\n')
        indent+='  '
        #GLOBAL
        out.write(indent + '<global>\n')
        indent+='  '
        out.write(indent + '<nt>64</nt>\n')
        out.write(indent + '<diskVectorDir>' + TmpDir + '</diskVectorDir>\n')
        out.write(indent + '<output>' + OutputDir + '</output>\n')
        out.write(indent + '<trajCounter>\n')
        indent+='  '
        out.write(indent + '<start>' + str(conf) + '</start>\n')
        out.write(indent + '<end>' + str(conf + 1) + '</end>\n')
        out.write(indent + '<step>40</step>\n')
        indent=indent[:-2]
        out.write(indent + '</trajCounter>\n')
        indent=indent[:-2]
        out.write(indent + '</global>\n')
        #/GLOBAL
        #A2AMATRIX
        out.write(indent + '<a2aMatrix>\n')
        indent+='  '
        #LOCAL FIELDS
        Nt=64
        interlace=4
        offset = 0
        nQuarks=1
        nDeltaT=4
        deltaT=[12,16,20,24]
        nG=4
        G=['GammaXGamma5','GammaYGamma5','GammaZGamma5','GammaTGamma5']
        Gname=['gx','gy','gz','g0']

        #nDeltaT=1
        #interlace=32
        #nG=1
        for i in range(0,nQuarks):
            for j in range(0,nQuarks):
                for t in range(0,Nt/interlace):
                    for dt in range(0,nDeltaT):
                        for ig in range(0,nG):
                            for ip in range(0,nMom):
                                tsrc=str(t*interlace+offset)
                                tsrc2=str((t*interlace+offset+deltaT[dt])%Nt)
                                out.write(indent + '<elem>\n')
                                indent+='  '
                                out.write(indent + '<file>' + PrefixUnsmeared + quarkL[i]+ '_' + str(tsrc) +'_'+quarkR[j]+'_' + str(tsrc2)+'.@traj@/'+G[ig]+'_'+moms[ip]+'.h5</file>\n')
                                out.write(indent + '<dataset>'+G[ig]+'_'+moms[ip]+'</dataset>\n')
                                out.write(indent + '<cacheSize>1</cacheSize>\n')
                                out.write(indent + '<name>'+Gname[ig]+'_'+str(deltaT[dt])+'_'+quarkL[i]+ '_'+quarkR[j]+'_' + str(tsrc) + '_p_' + str(ip) +'</name>\n')
                                indent=indent[:-2]
                                out.write(indent + '</elem>\n')
        #SINKS - THESE STAY AT [000]
        for i in range(0,nQuarks):
            for t in range(0,Nt/interlace):
                for dt in range(0,nDeltaT):
                    tsrc=str(t*interlace+offset)
                    tsrc2=str((t*interlace+offset+deltaT[dt])%Nt)
                    out.write(indent + '<elem>\n')
                    indent+='  '
                    out.write(indent + '<file>' + DirectorySmeared + 'Rho_'+str(tsrc2)+'_Phi_'+quark3[i]+ '_' + str(tsrc) +'.@traj@/Gamma5_0_0_0.h5</file>\n')
                    out.write(indent + '<dataset>Gamma5_0_0_0</dataset>\n')
                    out.write(indent + '<cacheSize>1</cacheSize>\n')
                    out.write(indent + '<name>sink_0_'+str(deltaT[dt])+'_'+quark3[i]+'_'+ str(tsrc) +'</name>\n')
                    indent=indent[:-2]
                    out.write(indent + '</elem>\n')
        for i in range(0,nQuarks):
            for t in range(0,Nt/interlace):
                for dt in range(0,nDeltaT):
                    for ip in range(0,nMom):
                        tsrc=str(t*interlace+offset)
                        tsrc2=str((t*interlace+offset+deltaT[dt])%Nt)
                        out.write(indent + '<elem>\n')
                        indent+='  '
                        out.write(indent + '<file>' + DirectorySmeared + 'Rho_'+str(tsrc2)+'_Phi_'+quark3[i]+ '_' + str(tsrc) +'.@traj@/Gamma5_'+moms_neg[ip]+'.h5</file>\n')
                        out.write(indent + '<dataset>Gamma5_'+moms_neg[ip]+'</dataset>\n')
                        out.write(indent + '<cacheSize>1</cacheSize>\n')
                        out.write(indent + '<name>sink_p_' + str(ip) + '_'+str(deltaT[dt])+'_'+quark3[i]+'_'+ str(tsrc) +'</name>\n')
                        indent=indent[:-2]
                        out.write(indent + '</elem>\n')
        #SINKS - THESE STAY AT [000] - AX
        for i in range(0,nQuarks):
            for t in range(0,Nt/interlace):
                for dt in range(0,nDeltaT):
                    tsrc=str(t*interlace+offset)
                    tsrc2=str((t*interlace+offset+deltaT[dt])%Nt)
                    out.write(indent + '<elem>\n')
                    indent+='  '
                    out.write(indent + '<file>' + DirectorySmeared + 'Rho_'+str(tsrc2)+'_Phi_'+quark3[i]+ '_' + str(tsrc) +'.@traj@/GammaTGamma5_0_0_0.h5</file>\n')
                    out.write(indent + '<dataset>GammaTGamma5_0_0_0</dataset>\n')
                    out.write(indent + '<cacheSize>1</cacheSize>\n')
                    out.write(indent + '<name>sink_ax_0_'+str(deltaT[dt])+'_'+quark3[i]+'_'+ str(tsrc) +'</name>\n')
                    indent=indent[:-2]
                    out.write(indent + '</elem>\n')
        for i in range(0,nQuarks):
            for t in range(0,Nt/interlace):
                for dt in range(0,nDeltaT):
                    for ip in range(0,nMom):
                        tsrc=str(t*interlace+offset)
                        tsrc2=str((t*interlace+offset+deltaT[dt])%Nt)
                        out.write(indent + '<elem>\n')
                        indent+='  '
                        out.write(indent + '<file>' + DirectorySmeared + 'Rho_'+str(tsrc2)+'_Phi_'+quark3[i]+ '_' + str(tsrc) +'.@traj@/GammaTGamma5_'+moms_neg[ip]+'.h5</file>\n')
                        out.write(indent + '<dataset>GammaTGamma5_'+moms_neg[ip]+'</dataset>\n')
                        out.write(indent + '<cacheSize>1</cacheSize>\n')
                        out.write(indent + '<name>sink_ax_p_' + str(ip) + '_'+str(deltaT[dt])+'_'+quark3[i]+'_'+ str(tsrc) +'</name>\n')
                        indent=indent[:-2]
                        out.write(indent + '</elem>\n')

        #2PT SINKS
        for i in range(0,nQuarks):
            for t in range(0,Nt/interlace):
                for ip in range(0,nMom):
                    tsrc=str(t*interlace+offset)
                    out.write(indent + '<elem>\n')
                    indent+='  '
                    out.write(indent + '<file>' + DirectorySmeared + 'Phi_'+quarkL[i]+ '_'+ str(tsrc)+'_Phi_'+quark3[i]+ '_' + str(tsrc) +'.@traj@/Identity_'+moms[ip]+'.h5</file>\n')
                    out.write(indent + '<dataset>Identity_'+moms[ip]+'</dataset>\n')
                    out.write(indent + '<cacheSize>1</cacheSize>\n')
                    out.write(indent + '<name>sink_p_' +str(ip) + '_'+quark3[i]+ '_'+ str(tsrc) +'</name>\n')
                    indent=indent[:-2]
                    out.write(indent + '</elem>\n')
        for i in range(0,nQuarks):
            for t in range(0,Nt/interlace):
                tsrc=str(t*interlace+offset)
                out.write(indent + '<elem>\n')
                indent+='  '
                out.write(indent + '<file>' + DirectorySmeared + 'Phi_'+quarkL[i]+ '_'+ str(tsrc)+'_Phi_'+quark3[i]+ '_' + str(tsrc) +'.@traj@/Identity_0_0_0.h5</file>\n')
                out.write(indent + '<dataset>Identity_0_0_0</dataset>\n')
                out.write(indent + '<cacheSize>1</cacheSize>\n')
                out.write(indent + '<name>sink_0_'+quark3[i]+ '_'+ str(tsrc) +'</name>\n')
                indent=indent[:-2]
                out.write(indent + '</elem>\n')
        #2PT SINKS - ax
        for i in range(0,nQuarks):
            for t in range(0,Nt/interlace):
                for ip in range(0,nMom):
                    tsrc=str(t*interlace+offset)
                    out.write(indent + '<elem>\n')
                    indent+='  '
                    out.write(indent + '<file>' + DirectorySmeared + 'Phi_'+quarkL[i]+ '_'+ str(tsrc)+'_Phi_'+quark3[i]+ '_' + str(tsrc) +'.@traj@/GammaT_'+moms[ip]+'.h5</file>\n')
                    out.write(indent + '<dataset>GammaT_'+moms[ip]+'</dataset>\n')
                    out.write(indent + '<cacheSize>1</cacheSize>\n')
                    out.write(indent + '<name>sink_ax_p_' +str(ip) + '_'+quark3[i]+ '_'+ str(tsrc) +'</name>\n')
                    indent=indent[:-2]
                    out.write(indent + '</elem>\n')
        for i in range(0,nQuarks):
            for t in range(0,Nt/interlace):
                tsrc=str(t*interlace+offset)
                out.write(indent + '<elem>\n')
                indent+='  '
                out.write(indent + '<file>' + DirectorySmeared + 'Phi_'+quarkL[i]+ '_'+ str(tsrc)+'_Phi_'+quark3[i]+ '_' + str(tsrc) +'.@traj@/GammaT_0_0_0.h5</file>\n')
                out.write(indent + '<dataset>GammaT_0_0_0</dataset>\n')
                out.write(indent + '<cacheSize>1</cacheSize>\n')
                out.write(indent + '<name>sink_ax_0_'+quark3[i]+ '_'+ str(tsrc) +'</name>\n')
                indent=indent[:-2]
                out.write(indent + '</elem>\n')
        #SOURCES
        for t in range(0,Nt/interlace):
            for ip in range(0,nMom):
                tsrc=str(t*interlace+offset)
                out.write(indent + '<elem>\n')
                indent+='  '
                out.write(indent + '<file>' + DirectorySmeared + 'Rho_'+ str(tsrc)+'_Rho_' + str(tsrc) +'.@traj@/Identity_'+moms_neg[ip]+'.h5</file>\n')
                out.write(indent + '<dataset>Identity_'+moms_neg[ip]+'</dataset>\n')
                out.write(indent + '<cacheSize>1</cacheSize>\n')
                out.write(indent + '<name>source_p_' + str(ip) + '_'+ str(tsrc) +'</name>\n')
                indent=indent[:-2]
                out.write(indent + '</elem>\n')
        for t in range(0,Nt/interlace):
            tsrc=str(t*interlace+offset)
            out.write(indent + '<elem>\n')
            indent+='  '
            out.write(indent + '<file>' + DirectorySmeared + 'Rho_'+ str(tsrc)+'_Rho_' + str(tsrc) +'.@traj@/Identity_0_0_0.h5</file>\n')
            out.write(indent + '<dataset>Identity_0_0_0</dataset>\n')
            out.write(indent + '<cacheSize>1</cacheSize>\n')
            out.write(indent + '<name>source_0_'+ str(tsrc) +'</name>\n')
            indent=indent[:-2]
            out.write(indent + '</elem>\n')
        #SOURCES - AX
        for t in range(0,Nt/interlace):
            for ip in range(0,nMom):
                tsrc=str(t*interlace+offset)
                out.write(indent + '<elem>\n')
                indent+='  '
                out.write(indent + '<file>' + DirectorySmeared + 'Rho_'+ str(tsrc)+'_Rho_' + str(tsrc) +'.@traj@/GammaT_'+moms_neg[ip]+'.h5</file>\n')
                out.write(indent + '<dataset>GammaT_'+moms_neg[ip]+'</dataset>\n')
                out.write(indent + '<cacheSize>1</cacheSize>\n')
                out.write(indent + '<name>source_ax_p_' + str(ip) + '_'+ str(tsrc) +'</name>\n')
                indent=indent[:-2]
                out.write(indent + '</elem>\n')
        for t in range(0,Nt/interlace):
            tsrc=str(t*interlace+offset)
            out.write(indent + '<elem>\n')
            indent+='  '
            out.write(indent + '<file>' + DirectorySmeared + 'Rho_'+ str(tsrc)+'_Rho_' + str(tsrc) +'.@traj@/GammaT_0_0_0.h5</file>\n')
            out.write(indent + '<dataset>GammaT_0_0_0</dataset>\n')
            out.write(indent + '<cacheSize>1</cacheSize>\n')
            out.write(indent + '<name>source_ax_0_'+ str(tsrc) +'</name>\n')
            indent=indent[:-2]
            out.write(indent + '</elem>\n')
        indent=indent[:-2]
        out.write(indent + '</a2aMatrix>\n')
        #/A2AMATRIX
        #PRODUCT
        out.write(indent + '<product>\n')
        indent+='  '
        #3PT
        if p2 != 2 or ContractFor != 2:
            for i in range(0,nQuarks):
                for j in range(0,nQuarks):
                    for t in range(0,Nt/interlace):
                        for dt in range(0,nDeltaT):
                            for ig in range(0,nG):
                                for ip in range(0,nMom):
                                    tsrc=str(t*interlace+offset)
                                    tsrc2=str((t*interlace+offset+deltaT[dt])%Nt)
                                    out.write(indent + '<elem>\n')
                                    indent+='  '
                                    out.write(indent + '<terms>sink_p_' + str(ip) + '_'+str(deltaT[dt])+'_'+quark3[i]+'_'+ str(tsrc) +' source_0_'+ str(tsrc) +' '+Gname[ig]+'_'+str(deltaT[dt])+'_'+quarkL[i]+ '_'+quarkR[j]+'_' + str(tsrc) +'_p_'+str(ip)+'</terms>\n')
                                    out.write(indent + '<times>\n')
                                    indent+='  '
                                    out.write(indent + '<elem>'+ str(tsrc2) +'</elem>\n')
                                    out.write(indent + '<elem>'+ str(tsrc) +'</elem>\n')
                                    indent=indent[:-2]
                                    out.write(indent + '</times>\n')
                                    out.write(indent + '<translations>0..63</translations>\n')
                                    out.write(indent + '<translationAverage>true</translationAverage>\n')
                                    indent=indent[:-2]
                                    out.write(indent + '</elem>\n')
        nAx=2
        axDescr=['','_ax']
        for iAx in range(0,nAx):
            for jAx in range(0,nAx):
                #if (iAx == 0 and jAx==0): continue
                #if (iAx == 0 or jAx==1): continue
                if p2 != 2 or ContractFor == 0:
                    #3PT - 0P
                    for i in range(0,nQuarks):
                        for j in range(0,nQuarks):
                            for t in range(0,Nt/interlace):
                                for dt in range(0,nDeltaT):
                                    for ig in range(0,nG):
                                        for ip in range(0,nMom):
                                            tsrc=str(t*interlace+offset)
                                            tsrc2=str((t*interlace+offset+deltaT[dt])%Nt)
                                            out.write(indent + '<elem>\n')
                                            indent+='  '
                                            out.write(indent + '<terms>sink'+axDescr[iAx]+'_0_'+str(deltaT[dt])+'_'+quark3[i]+'_'+ str(tsrc) +' source'+axDescr[jAx]+'_p_' + str(ip) + '_'+ str(tsrc) +' '+Gname[ig]+'_'+str(deltaT[dt])+'_'+quarkL[i]+ '_'+quarkR[j]+'_' + str(tsrc) +'_p_'+str(ip)+'</terms>\n')
                                            out.write(indent + '<times>\n')
                                            indent+='  '
                                            out.write(indent + '<elem>'+ str(tsrc2) +'</elem>\n')
                                            out.write(indent + '<elem>'+ str(tsrc) +'</elem>\n')
                                            indent=indent[:-2]
                                            out.write(indent + '</times>\n')
                                            out.write(indent + '<translations>0..63</translations>\n')
                                            out.write(indent + '<translationAverage>true</translationAverage>\n')
                                            indent=indent[:-2]
                                            out.write(indent + '</elem>\n')
                if p2 != 2 or ContractFor == 1:
                    #3PT - P0
                    for i in range(0,nQuarks):
                        for j in range(0,nQuarks):
                            for t in range(0,Nt/interlace):
                                for dt in range(0,nDeltaT):
                                    for ig in range(0,nG):
                                        for ip in range(0,nMom):
                                            tsrc=str(t*interlace+offset)
                                            tsrc2=str((t*interlace+offset+deltaT[dt])%Nt)
                                            out.write(indent + '<elem>\n')
                                            indent+='  '
                                            out.write(indent + '<terms>sink'+axDescr[iAx]+'_p_' + str(ip) + '_'+str(deltaT[dt])+'_'+quark3[i]+'_'+ str(tsrc) +' source'+axDescr[jAx]+'_0_'+ str(tsrc) +' '+Gname[ig]+'_'+str(deltaT[dt])+'_'+quarkL[i]+ '_'+quarkR[j]+'_' + str(tsrc) +'_p_'+str(ip)+'</terms>\n')
                                            out.write(indent + '<times>\n')
                                            indent+='  '
                                            out.write(indent + '<elem>'+ str(tsrc2) +'</elem>\n')
                                            out.write(indent + '<elem>'+ str(tsrc) +'</elem>\n')
                                            indent=indent[:-2]
                                            out.write(indent + '</times>\n')
                                            out.write(indent + '<translations>0..63</translations>\n')
                                            out.write(indent + '<translationAverage>true</translationAverage>\n')
                                            indent=indent[:-2]
                                            out.write(indent + '</elem>\n')
                if p2 != 2 or ContractFor == 2:
                    #2PT
                    for i in range(0,nQuarks):
                        for t in range(0,Nt/interlace):
                            for ip in range(0,nMom):
                                tsrc=str(t*interlace+offset)
                                out.write(indent + '<elem>\n')
                                indent+='  '
                                out.write(indent + '<terms>sink'+axDescr[iAx]+'_p_' + str(ip) + '_'+quark3[i]+ '_'+ str(tsrc) +' source'+axDescr[jAx]+'_p_' + str(ip) + '_'+ str(tsrc) +'</terms>\n')
                                out.write(indent + '<times>\n')
                                indent+='  '
                                out.write(indent + '<elem>'+ str(tsrc) +'</elem>\n')
                                indent=indent[:-2]
                                out.write(indent + '</times>\n')
                                out.write(indent + '<translations>0..63</translations>\n')
                                out.write(indent + '<translationAverage>true</translationAverage>\n')
                                indent=indent[:-2]
                                out.write(indent + '</elem>\n')
                    for i in range(0,nQuarks):
                        for t in range(0,Nt/interlace):
                            tsrc=str(t*interlace+offset)
                            out.write(indent + '<elem>\n')
                            indent+='  '
                            out.write(indent + '<terms>sink'+axDescr[iAx]+'_0_'+quark3[i]+ '_'+ str(tsrc) +' source'+axDescr[jAx]+'_0_'+ str(tsrc) +'</terms>\n')
                            out.write(indent + '<times>\n')
                            indent+='  '
                            out.write(indent + '<elem>'+ str(tsrc) +'</elem>\n')
                            indent=indent[:-2]
                            out.write(indent + '</times>\n')
                            out.write(indent + '<translations>0..63</translations>\n')
                            out.write(indent + '<translationAverage>true</translationAverage>\n')
                            indent=indent[:-2]
                            out.write(indent + '</elem>\n')
        indent=indent[:-2]
        out.write(indent + '</product>\n')
        #/PRODUCT
        indent=indent[:-2]
        out.write(indent + '</grid>\n')
        out.close()

        print 'xml created'
