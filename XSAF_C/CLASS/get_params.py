import numpy as np
import re
from pathlib import Path
from logging import getLogger

from regex import redshift_regex, redshift_rm_regex, transfer_file_rm_regex, regexs

lg=getLogger(__name__)

def generate_paramas_from_template(
    input_tempalte, output_tempalte_path, zrange, z, init_params
):
    if type(input_tempalte) == list:
        template = "\n".join(input_tempalte)
    else:
        with open(input_tempalte, "r") as file:
            template = file.read()

    # Initialisation of Camb parameters values for the linear case.
    template, edits = re.subn(redshift_rm_regex, "", template)
    template, edits = re.subn(transfer_file_rm_regex, "", template)
    for param, regex in regexs.items():
        template, edits = re.subn(regex, fr"\1 {init_params[param]}", template)
    redshift = ""
    for b in range(1, len(zrange) + 1):
        redshift += f"\ntransfer_redshift({b})    =  {z[len(zrange) - b]}"
    redshift += f"\ntransfer_redshift({len(zrange) + 1})    =  0.0"
    for b in range(1, len(zrange) + 2):
        redshift += f"\ntransfer_filename({b})    =  transfer_out{b}.dat"

    template, edits = re.subn(redshift_regex, fr"\1{redshift}", template)
    with open(output_tempalte_path, "w") as file:
        file.write(template)


def load_matterpower_data(root_prefix, zrange, is_class=False):
    use_cols = (
        0,
        1,
    )
    l=len(zrange)
    if is_class: # Get last params number (params00, params01, ...)
        prefixes={f.stem.split('_')[0] for f in Path('output').glob('params*')}
        last=max(prefixes) # Lexicographically
        lg.info(f'Last params: {last}')
    def fn(z):
        if is_class: # 1..398
            return f'output/{last}_z{z}_pk.dat'
        else:
            f"{root_prefix}_matterpower_{z}.dat"
        
    kh, pk = np.loadtxt( fn(l), usecols=use_cols, unpack=True )
    for a in range(l - 1, 0, -1):
        kh1, pk1 = np.loadtxt( fn(a), usecols=use_cols, unpack=True )
        kh = np.vstack((kh, kh1))
        pk = np.vstack((pk, pk1))
    kh1, pk1 = np.loadtxt( fn(l+1), usecols=use_cols, unpack=True )
    kh = np.vstack((kh, kh1))
    pk = np.vstack((pk, pk1))

    return kh, pk

def output(directory, zrange, z, pk, kh, out_params):
    with open(directory/"Pks8sqRatio_ist_LogSplineInterpPk.dat", "w") as outP:
        outP.write(f"{str('%.12e' % z[0])} to {str('%.12e' % z[len(z) - 1])}\n")
        outP.write(
            "{0[wb_new]} {0[h_new]} {0[wm_new]} {0[ns_new]} {0[wde_new]} {0[w0_new]} {0[wa_new]} {0[As_modified]}\n".format(out_params)
        )
        outP.write("   k(h/Mpc)      Pk/s8^2(Mpc/h)^3 \n")
        for files in range(len(zrange)):
            for flines in range(len(kh[files])):
                outP.write(
                    f"{str('%.12e' % kh[files][flines])} {str('%.12e' % pk[files][flines])}\n"
                )
