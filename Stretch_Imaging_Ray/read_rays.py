# -*- coding: utf-8 -*-
"""
Created on Thu May 26 11:44:52 2016

@author: cave
"""
import binascii
from struct import *
import numpy 
from numpy import frombuffer
import matplotlib.pyplot as plt





inBuffer = []


f = open("shift_shape_test_6_10.bin", "rb")

#try:
#    byte = f.read(8)
#    print byte
#    while byte != "":
#        # Do stuff with byte.
#        byte = f.read(8)
#        byte = bin(int(binascii.hexlify(byte), 16))
#        byte = byte[2:]
#        byte = int(byte, 2)
#        inBuffer.append(byte)
#        byte = f.read(8)
#finally:
#        
#    f.close()
#    
#Vddraw = inBuffer[3]*256|inBuffer[2]

try:
    #read perable data
    byte = f.read(4)
    byte_1 = numpy.frombuffer(byte, dtype=numpy.int32)
    number_of_rays=byte_1[0]
    print number_of_rays
    byte = f.read(4)
    byte_1 = numpy.frombuffer(byte, dtype=numpy.int32)
    lens_index=byte_1[0]
    byte = f.read(8)
    
    byte_1 = numpy.frombuffer(byte, dtype=numpy.float64)
    delta_at_len_index=byte_1[0]
    byte = f.read(8)
    byte_1 = numpy.frombuffer(byte, dtype=numpy.float64)
    lens_radius=byte_1[0]
    byte = f.read(8)
    byte_1 = numpy.frombuffer(byte, dtype=numpy.float64)
    curved_surface_y=byte_1[0]
    byte = f.read(8)
    byte_1 = numpy.frombuffer(byte, dtype=numpy.float64)
    sensor_width=byte_1[0]
    byte = f.read(8)
    byte_1 = numpy.frombuffer(byte, dtype=numpy.float64)
    aperture_width=byte_1[0]
    byte = f.read(8)
    byte_1 = numpy.frombuffer(byte, dtype=numpy.float64)
    aperture_height=byte_1[0]
    byte = f.read(8)
    byte_1 = numpy.frombuffer(byte, dtype=numpy.float64)
    delta_u=byte_1[0]
    byte = f.read(8)
    byte_1 = numpy.frombuffer(byte, dtype=numpy.float64)
    delta_v=byte_1[0]
    #preamble data read
    
    
    rays = np.zeros((number_of_rays, 5)) #Set up to read rays
    
    current_ray = 0 
    print number_of_rays
    while byte != "" and current_ray < number_of_rays:
        rays[current_ray,0] = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0] #Point 1 at the sensor
        rays[current_ray,1] = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0]
        rays[current_ray,2] = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0] #Point two, base of lens
        rays[current_ray,3] = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0]
        rays[current_ray,4] = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0] #Point three, curved surface of the lens
        if False: 
            print current_ray
            counter = (float(current_ray)/float(number_of_rays))*100
            print str(round(counter,2)) + "% Done"
        current_ray = current_ray + 1
        
finally:
    f.close()

    
