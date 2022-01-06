# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 10:49:37 2022

@author: marius bartkowski
"""

import pvlib
import os
import pandas as pd

df_IAM_mod_dir = []

for aoi in range(91):
    
    IAM_mod_dir = pvlib.iam.martin_ruiz(
            aoi,                                                       # The angle of incidence between the module normal vector and the sun-beam vector in degrees.
            a_r=0.25) 
    
    df_IAM_mod_dir.append(IAM_mod_dir)
    
df_IAM_mod_dir = pd.DataFrame(data=df_IAM_mod_dir)
    
df_IAM_mod_dir.to_excel(os.path.join("Exports", r'IAM.xlsx'))

