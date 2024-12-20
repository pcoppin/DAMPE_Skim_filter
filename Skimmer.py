# Python3 compatibility
from __future__ import print_function, division, absolute_import
import numpy as np

# Load DmpSoftware and ROOT
import ROOT
ROOT.gSystem.Load('libDmpEvent.so')
ROOT.gSystem.Load('libDmpService.so')

def TrueContainment(dc):
    nevents = dc.GetEntries()
    keys = ['pv_x', 'pv_y', 'pv_z', 'pvpart_px', 'pvpart_py', 'pvpart_pz']
    dd = {key: np.zeros(nevents, dtype=float) for key in keys}
    for i in range(nevents):
        ev = dc.GetDmpEvent(i)
        MC_truth = ev.pEvtSimuPrimaries()
        for key in keys:
            dd[key][i] = getattr(MC_truth, key)
    
    TopZ, BottomZ = -325, 448.
    cutTop, cutBottom = 440, 280
    topX = dd['pvpart_px']/dd['pvpart_pz']*(TopZ-dd['pv_z']) + dd['pv_x']
    topY = dd['pvpart_py']/dd['pvpart_pz']*(TopZ-dd['pv_z']) + dd['pv_y']
    bottomX = dd['pvpart_px']/dd['pvpart_pz']*(BottomZ-dd['pv_z']) + dd['pv_x']
    bottomY = dd['pvpart_py']/dd['pvpart_pz']*(BottomZ-dd['pv_z']) + dd['pv_y']
    w_fid = (abs(topX)<cutTop) * (abs(topY)<cutTop) * (abs(bottomX)<cutBottom) * (abs(bottomY)<cutBottom)
    ### REQUIRE THEM ALSO TO GO THROUGH THE TOP LAYER OF BGO!!!
    BottomZ = 44.
    cutBottom = 280
    bottomX = dd['pvpart_px']/dd['pvpart_pz']*(BottomZ-dd['pv_z']) + dd['pv_x']
    bottomY = dd['pvpart_py']/dd['pvpart_pz']*(BottomZ-dd['pv_z']) + dd['pv_y']
    w_fid = w_fid * (abs(bottomX)<cutBottom) * (abs(bottomY)<cutBottom)
    return w_fid

def SkimmerEv(dc, E_BGO_min=1):
    nevents = dc.GetEntries()
    skim = np.zeros(nevents, dtype=np.bool)
    for i in range(nevents):
        ev = dc.GetDmpEvent(i)

        # 1: BGO reconstructed energy greater than (15) GeV
        # Rationale for 15, we want to be able to cut on the quenching corrected energy at 20 GeV later
        # This way we reduce stastics without lossing necessary events
        BGOrec = ev.pEvtBgoRec()
        BGO_TotalE = BGOrec.GetTotalEnergy() # [MeV]
        if( BGO_TotalE<E_BGO_min ):
            continue

        # 2: rMaxELayer cut
        BGO_MaxELayer = max( [BGOrec.GetELayer(j) for j in range(14)] )
        ratio = BGO_MaxELayer/BGO_TotalE
        if( ratio>0.35 ):
            continue

        # 3: iBarMaxE cut
        #   Events where the maximum energy OF A SINGLE HIT is deposited at the edge of the calorimeter
        #   We skip the first layer of BGO as --> David: These events would prob be removed by vertex cut anyway
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
                ### The following 3 definitions are all equivalent!
                IdxBar = ROOT.DmpBgoBase.GetBarID(int(BGO_HitBarId[j]))
                # IdxBar = np.divide(BGO_HitBarId[j],2**6)%32
                # IdxBar = (BGO_HitBarId[j] >> 6) & 0x1f
                Hit123_maxE[HitLayer-1] = HitE
                Hit123_BarIdx[HitLayer-1] = IdxBar
        Criteria_is_met = False
        for j in range(3):
            if( Hit123_BarIdx[j]<=0 or Hit123_BarIdx[j]==21 ):
                Criteria_is_met = True
        if( Criteria_is_met ):
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