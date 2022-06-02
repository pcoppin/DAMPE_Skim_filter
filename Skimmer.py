# Python3 compatibility
from __future__ import print_function, division, absolute_import
import numpy as np

# Load DmpSoftware and ROOT
from ROOT import gSystem
gSystem.Load('libDmpEvent.so')
from ROOT import DmpChain

def SkimmerEv(dc):
    nevents = dc.GetEntries()
    skim = np.zeros(nevents, dtype=np.bool)
    for i in range(nevents):
        ev = dc.GetDmpEvent(i)

        # 1: BGO reconstructed energy greater than 20 GeV
        BGOrec = ev.pEvtBgoRec()
        BGO_TotalE = BGOrec.GetTotalEnergy() # [MeV]
        if( BGO_TotalE<2e4 ):
            continue

        # 2: rMaxELayer cut
        BGO_MaxELayer = max( [BGOrec.GetELayer(j) for j in range(14)] )
        ratio = BGO_MaxELayer/BGO_TotalE
        if( ratio>0.35 ):
            continue

        # 3: iBarMaxE cut
        #   Events where the maximum energy OF A SINGLE HIT is deposited at the edge of the calorimeter?
        #   By 1, 2, 3    So how come we skip the first layer of BGO? --> David: Would prob be removed by vertex cut anyway
        BGOhits = ev.pEvtBgoHits()
        BGO_NHits = BGOhits.GetHittedBarNumber()
        BGO_HitE = np.array(BGOhits.fEnergy, dtype=np.float)
        BGO_HitBarId = np.array(BGOhits.fGlobalBarID, dtype=np.uint16)

        Hit123_maxE = np.zeros(3)
        Hit123_BarIdx = -1*np.ones(3)
        for j in range(BGO_NHits):
            HitE = BGO_HitE[j]
            HitLayer = BGOhits.GetLayerID(j)
            if( (HitLayer in [1,2,3]) and (HitE>Hit123_maxE[HitLayer-1]) ):
                IdxBar = np.divide(BGO_HitBarId[j],2**6)%31
                Hit123_maxE[HitLayer-1] = HitE
                Hit123_BarIdx[HitLayer-1] = IdxBar
        for j in range(3):
            if( Hit123_BarIdx[j]<=0 or Hit123_BarIdx[j]==21 ):
                continue

        # 4 & 5: FullBGO cut (+ shower vertex/direction successfully reconstructed)
        BGO_ZTop, BGO_ZBottom = 46, 448
        BGO_Slope = [BGOrec.GetSlopeYZ(), BGOrec.GetSlopeXZ()]
        BGO_Intercept = [BGOrec.GetInterceptYZ(), BGOrec.GetInterceptXZ()]
        if( (BGO_Slope[1]==0 and BGO_Intercept[1]==0) or (BGO_Slope[0]==0 and BGO_Intercept[0]==0) ):
            continue
        XTop = BGO_Slope[1]*BGO_ZTop + BGO_Intercept[1]
        YTop = BGO_Slope[0]*BGO_ZTop + BGO_Intercept[0]
        XBottom = BGO_Slope[1]*BGO_ZBottom + BGO_Intercept[1]
        YBottom = BGO_Slope[0]*BGO_ZBottom + BGO_Intercept[0]
        if( np.any(np.abs([XTop,YTop,XBottom,YBottom])>280) ):
            continue

        # If we were not stopped by now, set the skim variable to true
        skim[i] = True
    return skim