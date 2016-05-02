import numpy as np

class instrument(object):
    def __init__(self, name, camnum):
        self.name = name
        self.camname  = ['C'+str(i) for i in np.arange(camnum)]
        
        self.flatname = 'flat'
        self.biasname = 'bias'
        self.darkname = 'dark'
        self.objname  = 'img'
        
        self.ftype_post = {self.objname: 'o', self.flatname: 'f', 
            self.biasname: 'b', self.darkname: 'd'}
        
        self.split = False
        self.cam_bias   = [True]
        self.cam_dark   = [True]
        self.filters    = []
        self.filterscam = []
        

# RATIR specific initialization
ratir = instrument('ratir', 4)
ratir.filters = ['u', 'g', 'r', 'i', 'Z', 'Y', 'J', 'H']
ratir.filterscam = ['C0', 'C0', 'C0', 'C1', 'C2', 'C2', 'C3','C3']
ratir.cam_bias = [True, True, False, False]
ratir.cam_dark = [True, True, False, False]

# RATIR specific due to software gain when saving files
ratir.cam_satur = [ lambda SOFTGAIN: (2.**16/SOFTGAIN)-1, lambda SOFTGAIN: (2.**16/SOFTGAIN)-1, lambda SOFTGAIN: (36000./SOFTGAIN)-1, lambda SOFTGAIN: (36000./SOFTGAIN)-1 ] # saturation levels for each detector in DNs as a function of the SOFTGAIN keyword extracted from a frame's header

# Only needed for split cameras
ratir.split = True
ratir.CAM_SPLIT = [False, False, True, True]
# RATIR H2RG filter slices
ratir.C0_SLICE = np.s_[0:1023,125:1023]
ratir.C1_SLICE = np.s_[75:1000,15:1000]
ratir.Z_SLICE = np.s_[4:975,100:2000]
ratir.Y_SLICE = np.s_[1135:2043,240:2043]
ratir.J_SLICE = np.s_[50:1000,4:2000]
ratir.H_SLICE = np.s_[1200:2043,4:1940]
#ratir.H2RG_SLICES = [ ratir.Z_SLICE, ratir.J_SLICE, ratir.Y_SLICE, ratir.H_SLICE ] # same order as H2RG_FILTERS.  0+2 are C2, 1+3 are C3

ratir.slice = {'C0': ratir.C0_SLICE, 'C1': ratir.C1_SLICE, 'C2a':ratir.Z_SLICE, 'C2b':ratir.Y_SLICE, 'C3a':ratir.J_SLICE, 'C3b':ratir.H_SLICE}
ratir.SLICE_FILTERS = {'C2a': 'Z', 'C2b': 'Y', 'C3a': 'J', 'C3b': 'H'}
#ratir.SLICES = {'Z': ratir.Z_SLICE, 'J': ratir.J_SLICE, 'Y': ratir.Y_SLICE, 'H': ratir.H_SLICE, 'C0': ratir.C0_SLICE, 'C1': ratir.C1_SLICE}
#ratir.SPLIT_FILTERS = [ 'Z', 'J', 'Y', 'H' ] # RATIR NIR bands.  0+2 are C2, 1+3 are C3


"""
PROPOSALS = { 'NIRstandard': '0000', 'OPTstandard': '00001', 'cluster': '0002', 'galaxy': '0003', 'blank': '0004', 'pointing': '0005', 'bias': '0006', 'dark': '0007', 'flat': '0008', 'focus': '0009', 'misc': '0010', 'GRB': '1000' } # 2012 proposal names and id numbers. source: rsync://ratir.astroscu.unam.mx/public/proposalidentifiers.txt
CAM_PXSCALE = [0.32, 0.32, 0.3, 0.3] # C0, C1, C2, C3 in arcsec/px


CONFIG_LOCATION = 'astro_functs.py' # name of file containing configuration information, currently this file.

CAM_WAVE  = ['OPT', 'OPT', 'IR', 'IR']
CAM_SPLIT = [False, False, True, True]
#WCS relevant parameters (RATIR H2RGs have barrel distortions)
a = -19.60381671
b = -4128.15179797
CAM_SECPIX1  = [0.3168, 0.3171, 0.2988, 0.2983]
CAM_SECPIX2  = [0.3171, 0.3191, 0.2955, -0.2945]
CAM_THETA    = [0.60, 2.40, -88.1, 90.4]
CAM_X0       = [512, 512, 1177, 924]
CAM_Y0       = [512, 512, 1031, 982]

H2RG_ASTR = {'PV1_1': 1.0, 'PV2_1': 1.0, 'PV1_17':a, 'PV2_17':a, 'PV1_19': 2.0*a, 'PV2_19':2.0*a, 'PV1_21':a, 'PV2_21':a, 'PV1_31':b, 'PV2_31': b, 'PV1_33':3.0*b, 'PV2_33':3.0*b, 'PV1_35':3.0*b, 'PV2_35':3.0*b, 'PV1_37':b, 'PV2_37': b}
RA_KEY = 'ETRRQRA'
DEC_KEY = 'ETRRQDE'
OFFRA_KEY = 'ETRRQRAO'
OFFDEC_KEY = 'ETRRQDEO'
SOFTGAIN_KEY = 'SOFTGAIN'


CAM_GAIN = [ lambda SOFTGAIN: 16.80/SOFTGAIN, lambda SOFTGAIN: 18.64/SOFTGAIN, lambda SOFTGAIN: 2.2/SOFTGAIN, lambda SOFTGAIN: 2.4/SOFTGAIN ] # gain of each camera as a function of the SOFTGAIN keyword extracted from a frame's header
CENTER_KEY = 'STRRQAP' # RATIR header keyword specifying which H2RG filters the target is focused on
# frame corners in arcmin offset from center.  top-left, bottom-left, bottom-right, top-right.  top==north, left==east
CAMOFFS = np.array([[[2.785,2.632], [2.604,-2.775], [-2.800,-2.615], [-2.635,2.789]],   # C0 corner offsets in arcmin
                    [[2.817,2.624], [2.607,-2.818], [-2.807,-2.624], [-2.616,2.818]],   # C1 corner offsets in arcmin
                    [[5.012,6.229], [4.569,-3.905], [-0.227,-3.678], [0.228,6.453]],    # C2-Z corner offsets in arcmin
                    [[-0.556,6.488],[-1.013,-3.642],[-5.488,-3.445], [-5.018,6.683]],   # C2-Y corner offsets in arcmin
                    [[4.834,5.430], [4.720,-4.701], [-0.324,-4.647], [-0.185,5.494]],   # C3-J corner offsets in arcmin
                    [[-0.916,5.503],[-1.059,-4.639],[-5.318,-4.594], [-5.154,5.557]]])  # C3-H corner offsets in arcmin
FRAMECENTER = CAMOFFS.mean(axis=1) # field centers in arcmin (E,N) offset from center
# Offsets of the pointing apertures east and north in arcmin
APOFFS  = { "rcenter":      np.array([0,0]),
            "icenter":      np.array([0,0]),
            "ricenter":     np.array([0,0]),
            "riZJcenter":   np.array([1.2,0]),
            "riYHcenter":   np.array([-1.8,0]),
            "ZJcenter":     np.array([2.2,0.7]),
            "YHcenter":     np.array([-3.2,0.7])}
  
"""        
