#!/usr/bin/env python

###################################################################################################
# Authors: S. Yahia-Cherif
# Last update: 29/06/2020
# This script call Camb to get the matter power spectrums.
###################################################################################################

import os
import numpy as np
import pandas as pd
import sys, platform
from time import time
from get_params import load_matterpower_data, output
from integral import get_integral
from clear import clear
from pathlib import Path
import re
from logging import basicConfig, getLogger
from argparse import ArgumentParser
from dataclasses import dataclass, field
from typing import List
from json import loads
from copy import deepcopy

template_fn='explanatory_base.ini'
qtl_dir=Path('../../QTLauncher')
d_step = 3

ap=ArgumentParser(description='')
ap.add_argument('-l', dest='loglevel', help='Loglevel', default='INFO')
args=ap.parse_args()

basicConfig(format='{asctime}: {message}', datefmt='%d.%m %H:%M:%S', style='{',
            level=args.loglevel)
lg=getLogger(__name__)

clear(is_class=True)
s = time()

# Loading the parameters files.
Codes_elts = np.loadtxt(qtl_dir/"Codes_W.txt", dtype="str")
X_elts = np.loadtxt(qtl_dir/"XSAF_W.txt", dtype="str")
Parameters_elts = np.loadtxt(qtl_dir/"Parameters_W.txt", dtype="str")

# df=pd.read_table(qtl_dir/"Codes_W.txt", sep='\s+', header=None, index_col=0)
# loads(df.to_json(orient="columns"))['1']

# All the parameters are saved in a dictionary.
XSAF_elts = {}
for i in range(len(X_elts)):
    XSAF_elts[X_elts[i][0]] = X_elts[i][1]
for i in range(len(Parameters_elts)):
    XSAF_elts[Parameters_elts[i][0]] = Parameters_elts[i][1]

# Initializing all the variables.
fid = np.array(
    [
        float(XSAF_elts["VFidOmegab"]),
        float(XSAF_elts["VFidh"]),
        float(XSAF_elts["VFidOmegam"]),
        float(XSAF_elts["VFidns"]),
        float(XSAF_elts["VFidOmegaDE"]),
        float(XSAF_elts["VFidw0"]),
        float(XSAF_elts["VFidwa"]),
        float(XSAF_elts["VFidAs"]) * 10 ** (-9),
    ]
)

z = np.linspace(
    float(XSAF_elts["Vzmin"]), float(XSAF_elts["Vzmax"]), int(XSAF_elts["Vprec_Int_z"])
)
z_win_3_1 = 1.0 / 4 * (3 * z[1:] + z[:-1])
z_win_1_3 = 1.0 / 4 * (z[1:] + 3 * z[:-1])
z_win_1_2 = 1.0 / 2 * (z[1:] + z[:-1])
z_win = np.zeros(len(z) + len(z_win_3_1) + len(z_win_1_2) + len(z_win_1_3))
i, j = 0, 0
while i < len(z_win):
    if i == 0:
        z_win[i] = z[j]
        i = i + 1
        j = j + 1
    elif i > 0 and i < len(z_win):
        z_win[i] = z_win_1_3[j - 1]
        z_win[i + 1] = z_win_1_2[j - 1]
        z_win[i + 2] = z_win_3_1[j - 1]
        z_win[i + 3] = z[j]
        i = i + 4
        j = j + 1
z = np.copy(z_win)

zrange = range(len(z))
z_str = np.chararray((len(zrange)), itemsize=1000)

z = z[0 : len(zrange)]
zrange = zrange[0 : len(zrange)]

z_str = np.chararray((len(zrange)), itemsize=1000)
i = 0
while i < len(zrange):
    z_str[i] = "iz_" + str(zrange[i])
    i = i + 1
i = 0
z_str = z_str.decode("utf-8")

w_neu = float(XSAF_elts["VFidWnu"])

path_repository = Path("WP_Pk_seq_regex")

#     for _ in {'wb', 'h', 'wm', 'ns', 'wde' ,'w0', 'wa', 's8'}

@dataclass
class pair:
    val: float
    eps: float

    def __str__(self):
        return str(self.val)
    
@dataclass
class params:
    wb: pair
    h: pair
    wm: pair
    ns: pair
    wde: pair
    w0: pair
    wa: pair
    s8: pair

    # Not used, here we have omega_x and Omega_x.
    # def use_physical(self, f):
    #     '''omega to omega?'''
    #     return float(f) * self.h.val ** 2

# VFidOmegab 0.05
# VFidh 0.67
# VFidOmegam 0.32
# VFidw0 -1
# VFidwa 0
# VFidns 0.96
# VFidOmegaDE 0.68
# VFidsigma8 0.815534

    # .ini variable names
    def create_params_ini(self, do_nonlinear, scalar_amp, output_root="test"):
        params={
            "Omega_b": self.wb.val,
            "Omega_m": self.wm.val,
            "h": self.h.val,
            "n_s": self.ns.val,            
            # "omega_lambda": 1.0 - wm_new,
            "scalar_amp": scalar_amp,
            "output_root": output_root,
        }
        if XSAF_elts['UseOmegaDE_ch'] == '0':
            params.update({
                'fluid_equation_of_state': 'EDE',
                "Omega_EDE": self.wde.val
            })
        else:
            params.update({
                'fluid_equation_of_state': 'CLP',
                "w0_fld": self.w0.val,
                "wa_fld": self.wa.val
            })
        
        # Now convert to the .ini format
        if do_nonlinear:
            params['non_linear']='' # hwcode        
        else:
            pass # Leave empty
        
        # Now write params.ini
        lg.debug(f'Creating params.ini with {params=}')
        res_params={ **all_params, **params }
        with open('params.ini', 'w') as fo:
            for k, v in res_params.items(): fo.write(f'{k}={v}\n')

    def out_params(self, As_modified):
        return {
            "wb_new": self.wb,
            "h_new": self.h,
            "wm_new": self.wm,
            "ns_new": self.ns,
            "wde_new": self.wde,
            "w0_new": self.w0,
            "wa_new": self.wa,
            "As_modified": As_modified
        }
        
base_params=params(
    wb=pair(fid[0], float(XSAF_elts["VStepOmegab"])),
    h=pair(fid[1], float(XSAF_elts["VSteph"])),
    wm=pair(fid[2], float(XSAF_elts["VStepOmegam"])),
    ns=pair(fid[3], float(XSAF_elts["VStepns"])),
    wde=pair(fid[4], float(XSAF_elts["VStepOmegaDE"])),
    w0=pair(fid[5], float(XSAF_elts["VStepw0"])),
    wa=pair(fid[6], float(XSAF_elts["VStepwa"])),
    s8=pair(fid[7], float(XSAF_elts["VStepsigma8"]))
)

# The path are all declared and initialized in variables.
fold_path_fid=path_repository/"fid"

# The directories are all created if they don't exist.
fold_path_fid.mkdir(exist_ok=True)

# for k, v in folders.items():
#     for f in v: f.mkdir(exist_ok=True)

all_params={}
with open(template_fn) as fi:
    for line in fi.readlines():
        ls=line.strip()
        if not ls or ls[0] == '#': continue # Skip comments

        if 0: # There are parameters with spaces
            if m := re.match('^(\w+)\s*=\s*(.*)$', ls):
                all_params[gr[0]]=(gr:=m.groups())[1]
            else:
                breakpoint()
        else:
            lp=ls.partition('=')
            all_params[lp[0].strip()]=lp[2].strip()

    print(f'{len(all_params)=}')

def run_class():
    res=os.system("./class params.ini >output/class.log 2>&1")
    assert not res

# step
st=int((d_step - 1) / 2) # Doesn't change in the loop

suffixes={ # Depending on step
    -3: 'dw3',
    -2: 'dw2',
    -1: 'dw',
    1: 'up',
    2: 'up2',
    3: 'up3',
}
runs=('main', 'wb', 'h', 'wm', 'ns', 'wde' ,'w0', 'wa', 's8')
#runs=('main', 's8')
for p in runs:
    # Variations with different parameters, changing one at a time
    lg.info(f'Running for {p}')
    
    for l in range(-st, st+1): # Few steps to each side
        if not l: continue # No offset

        lg.info(f'{l=}')

        if p in {'main', 's8'}:
            params=base_params # Won't change it
        else:
            params=deepcopy(base_params)
            par=getattr(params, p)
            par.val+=l*par.eps

        # linear case (non fiducial).
        params.create_params_ini(0, params.s8.val)
        #os.system(f'cp params.ini {path_repository}/params_{p}{l:+}_l.ini')
        run_class()
        
        kh, pk = load_matterpower_data("test", zrange, is_class=True)
        # Interpolate and compute the integrals of the power spectrums.
        integrale = get_integral(pk, kh)

        # Compute sigma8 and As.
        if p == 'main':
            # Compute sigma8 and As. And then they'll be modified on each step
            s8_fiducial_ref = np.sqrt(integrale)
            As_fiducial_ref = np.copy(params.s8.val)

            As_modified=As_fiducial_ref # For the output
        elif p == 's8':
            s8_modified = s8_fiducial_ref * (1 + l * params.s8.eps)
            As_modified = As_fiducial_ref * (s8_modified / s8_fiducial_ref) ** 2
        else:
            s8_modified = np.sqrt(integrale)
            As_modified = As_fiducial_ref * (s8_fiducial_ref / s8_modified) ** 2

        # non linear case (non fiducial).
        lg.debug(f'{As_modified=}')
        params.create_params_ini(1, As_modified) # Will create similar params
        #os.system(f'cp params.ini {path_repository}/params_{p}{l:+}_nl.ini')
        run_class()
        
        kh, pk = load_matterpower_data("test", zrange, is_class=True)

        if p == 'main':
            kh = np.delete(kh, len(kh) - 1, axis=0)
            pk = np.delete(pk, len(pk) - 1, axis=0)

            dn=fold_path_fid
        else:
            # Save all the non fiducial power spectrums.
            dn=path_repository/f'{p}_{suffixes[l]}'
            dn.mkdir(exist_ok=True)

        #As=params.s8.val if p == 'main' else As_modified
        output(dn, zrange, z, pk, kh, params.out_params(As_modified))

        if p == 'main': break # No steps here

# Erase all the temporar files.
clear(is_class=True)
print("########## FINISH:", time() - s)
