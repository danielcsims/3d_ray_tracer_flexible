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
    byte = f.read(4)
    byte_1 = numpy.frombuffer(byte, dtype=numpy.int32)
    print byte_1
    height=byte_1[0]
    byte = f.read(4)
    byte_1 = numpy.frombuffer(byte, dtype=numpy.int32)
    print byte_1
    width=byte_1[0]
    recorded_array = np.zeros((height,width))
    current_column = 0
    current_row = 0
    while byte != "":
        # Do stuff with byte.
        byte = f.read(8)
        byte_1 = numpy.frombuffer(byte, dtype=numpy.float64)
        if byte != "":
            recorded_array[current_row-1,current_column-1]=byte_1
        if current_column != width-1:
            current_column=current_column+1
        else:
            current_column = 0
            current_row = current_row+1
        if byte_1 != 0.0:
            print byte_1
        if current_row == height:
            byte = ""
    result = recorded_array[149,:]
    byte = f.read(4)
    sqrt_of_number_of_points = numpy.frombuffer(byte, dtype=numpy.int32)[0]
    print sqrt_of_number_of_points
    byte = f.read(8)
    lens_radius = numpy.frombuffer(byte, dtype=numpy.float64)[0]
    print lens_radius
    byte = f.read(8)
    circle_height = numpy.frombuffer(byte, dtype=numpy.float64)[0]
    print circle_height
    current_point = 0
    while byte != "" and current_point < (sqrt_of_number_of_points*sqrt_of_number_of_points):
        x1 = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0] #Point 1 at the sensor
        y1 = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0]
        x2 = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0] #Point two, base of lens
        y2 = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0]
        x3 = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0] #Point three, curved surface of the lens
        y3 = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0]
        x4 = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0] #Point four, 
        y4 = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0]
        if True:
            #print "---"
            #print x1
            #print y1
            #print x2
            #print y2
            #print x3
            #print y3
            #print x4
            #print y4
            print "\draw (%f,%f) -- (%f,%f) -- (%f,%f) -- (%f,%f) [line width=0.1mm];" % (x1,y1,x2,y2,x3,y3,x4,y4)
        current_point = current_point + 1
        
finally:
    f.close()

    
