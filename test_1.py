# -*- coding: utf-8 -*-
"""
Created on Thu May 26 11:44:52 2016

@author: cave
"""

f = open("myfile.bin.bin", "rb")
try:
    byte = f.read(1)
    while byte != "":
        # Do stuff with byte.
        byte = f.read(1)
        print byte
finally:
        
    f.close()