import sys
import os
import re

pre_CC_path = ["Cl_GG/", "Cl_LL/", "Cl_GL/"]

# 2 modes : Standard with 10 bins ("S") / Extended with 11 bins ("E") / Other ("O") - Equidistant redshift (Higher FoM), Common Bias with 5 bins ("C")
#           since integration is better ("a rough explanation"))
model = 'E'
# 2 subcases : PESS = SEMI-PESS + zcut / PESS = SEMI_PESS without zcut / OPT
subCase = 'OPT'
# curvature : FLAT/NON_FLAT
curvature = 'F'
# Gamma : Yes or No
gamma = 'N'
# Zcut : Yes or No
zcut = 'N'
# Compute WP_Pk
Powspec_choice = 'Y'
# Compute Photo-z
photoZ = 'Y'
# Compute C_l's
compute_Cl = 'Y'
# Density photo: 30/10 = 3 gal.arcmin^-2 / bin for photo
#           0.35 gal.arcmin^2 on all bins (4 for IST and 0.35 for 5 here for spectro)
density = 30

# Store into file
with open('XCPP_params.txt', 'w') as f:
  f.write(model+'\n')
  f.write(subCase+'\n')
  f.write(curvature+'\n')
  f.write(gamma+'\n')
  f.write(zcut+'\n')
  f.write(photoZ+'\n')
  f.write(compute_Cl+'\n')

# Compute P_k
def Using_WP_Pk(Powspec_choice):
  if (Powspec_choice == 'Y'):
    print('Removing current WP_Pk CAMB')
    os.system('rm -rf ./WP_Pk')
    print('Calling CAMB')
    os.chdir('CAMB-0.1.7')
    os.system('rm -rf ./WP_Pk')
    os.system('python3 Camb_launcher_XSAF.py')
    os.system('cp -rf ./WP_Pk ../')
    os.chdir('../')

# Load params
if (model == 'S'):
  numBias = 10
  numParams = numBias + 12
elif (model == 'E'):
  numBias = 11
  numParams = numBias + 12
elif (model == 'O' or model =='C'):
  # CAUTION : to write as function of number of bins into params.txt
  if model == 'O':
    # Set any wanted value
    numBias = 10
    numParams = numBias + 12
  elif model == 'C':
    # Common bias
    numBias = 5
    numParams = numBias + 12
else: 
  print('Error on model O or C !')
  sys.exit()

# Update number of bins into params.txt and XSAF_W.txt
with open('params.txt', 'r') as f:
  data = f.read()
  data = re.sub('(VRS_bins).*\n', r'\1 '+str(numBias)+'\n',data)
  data = re.sub('(VGdensity).*\n', r'\1 '+str(density)+'\n',data)
  # CURVATURE
  if (curvature == 'F'):
    data = re.sub('(FNF_ch).*\n', r'\1 '+'0'+'\n',data)
    data = re.sub('(UseOmegaDE_ch).*\n', r'\1 '+'1'+'\n',data)
  elif (curvature == 'NF'):
    data = re.sub('(FNF_ch).*\n', r'\1 '+'1'+'\n',data)
    data = re.sub('(UseOmegaDE_ch).*\n', r'\1 '+'0'+'\n',data)
  else:
    print('Error on Curvature in params.txt')
  # NO-GAMMA / GAMMA
  if (gamma == 'N'):
    data = re.sub('(Usegamma_ch).*\n', r'\1 '+'1'+'\n',data)
  elif (gamma == 'Y'):
    data = re.sub('(Usegamma_ch).*\n', r'\1 '+'0'+'\n',data)
  else:
    print('Error on UseGamma in params.txt')

with open('params.txt', 'w') as f:
  f.write(data)
with open('../QTLauncher/XSAF_W.txt', 'w') as f:
  f.write(data)

# Computing WP_Pk
Using_WP_Pk(Powspec_choice)

# Compute Cl's
if (compute_Cl == 'Y'):
  os.system('rm -rf Cl_GG/ Cl_LL/ Cl_GL/')

for paramo in range(numParams):
	if paramo == 0:
	    CC_path = ["C_wb_up", "C_wb_up2", "C_fid", "C_wb_dw", "C_wb_dw2"]
	elif paramo == 1:
	    CC_path = ["C_h_up", "C_h_up2", "C_fid", "C_h_dw", "C_h_dw2"]
	elif paramo == 2:
	    CC_path = ["C_wm_up", "C_wm_up2", "C_fid", "C_wm_dw", "C_wm_dw2"]
	elif paramo == 3:
	    CC_path = ["C_ns_up", "C_ns_up2", "C_fid", "C_ns_dw", "C_ns_dw2"] 
	elif paramo == 4:
	    CC_path = ["C_wde_up", "C_wde_up2", "C_fid", "C_wde_dw", "C_wde_dw2"]  
	elif paramo == 5:
	    CC_path = ["C_w0_up", "C_w0_up2", "C_fid", "C_w0_dw", "C_w0_dw2"]
	elif paramo == 6:
	    CC_path = ["C_wa_up", "C_wa_up2", "C_fid", "C_wa_dw", "C_wa_dw2"]
	elif paramo == 7:
	    CC_path = ["C_s8_up", "C_s8_up2", "C_fid", "C_s8_dw", "C_s8_dw2"]
	elif paramo == 8:
	    CC_path = ["C_gamma_up", "C_gamma_up2", "C_fid", "C_gamma_dw", "C_gamma_dw2"]
	elif paramo == 9:
	    CC_path = ["C_A_IA_up", "C_A_IA_up2", "C_fid", "C_A_IA_dw", "C_A_IA_dw2"]
	elif paramo == 10:
	    CC_path = ["C_n_IA_up", "C_n_IA_up2", "C_fid", "C_n_IA_dw", "C_n_IA_dw2"]
	elif paramo == 11:
	    CC_path = ["C_B_IA_up", "C_B_IA_up2", "C_fid", "C_B_IA_dw", "C_B_IA_dw2"]
	elif paramo > 11:
	    CC_path = ["C_b"+str(paramo-11)+"_up", "C_b"+str(paramo-11)+"_up2", "C_fid", "C_b"+str(paramo-11)+"_dw", "C_b"+str(paramo-11)+"_dw2"]

	for i in range(len(CC_path)):
	    if not os.path.exists(pre_CC_path[0]+CC_path[i]):
	        os.makedirs(pre_CC_path[0]+CC_path[i])
	    if not os.path.exists(pre_CC_path[1]+CC_path[i]):
	        os.makedirs(pre_CC_path[1]+CC_path[i])
	    if not os.path.exists(pre_CC_path[2]+CC_path[i]):
	        os.makedirs(pre_CC_path[2]+CC_path[i])

# Create output directory
if not os.path.exists("output"):
   os.makedirs("output")
os.system('make gnu && printf "\\nCOMPILATION OK !\\n"')
os.system('/usr/bin/time -p ./main_gnu.exe')
