# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 12:22:08 2016

@author: sautt
"""

with open("ProjectParameters.json", 'r') as f:
    text = f.read()
    text = text.replace('"factor"', '"modulus"')
    text = text.replace('implemented_in_module', 'kratos_module')
    text = text.replace('implemented_in_file', 'python_module')
    f.close()

with open("ProjectParameters.json", 'w') as f:
    f.write(text)
    f.close()


with open("MainKratos.py", 'r') as k:
    text = k.read()
    text = text.replace('GetComputeModelPart()', 'GetComputingModelPart()')
    k.close()

with open("MainKratos.py", 'w') as k:
    k.write(text)
    k.close()


