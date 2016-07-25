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

f = open("test_debug_6_13_2_1000_try_2.bin", "rb")

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


PSF_resolution = 2048*2 #number of boxes that stores the PSF 
PSF_width = 60 #mm
distance_of_scene = 100 #mm 
try:
    byte = f.read(4)
    byte_1 = numpy.frombuffer(byte, dtype=numpy.int32)
    print byte_1
    sqrt_of_number_of_points=byte_1[0]
    byte = f.read(4)
    byte_1 = numpy.frombuffer(byte, dtype=numpy.int32)
    print byte_1
    lens_index=byte_1[0]
    byte = f.read(8)
    byte_1 = numpy.frombuffer(byte, dtype=numpy.double)
    print byte_1
    delta_at_len_index=byte_1[0]
    byte = f.read(8)
    byte_1 = numpy.frombuffer(byte, dtype=numpy.double)
    print byte_1
    lens_radius=byte_1[0]
    byte = f.read(8)
    byte_1 = numpy.frombuffer(byte, dtype=numpy.double)
    print byte_1
    curved_surface_y=byte_1[0]
    byte = f.read(8)
    byte_1 = numpy.frombuffer(byte, dtype=numpy.double)
    print byte_1
    sensor_width=byte_1[0]
    byte = f.read(8)
    byte_1 = numpy.frombuffer(byte, dtype=numpy.double)
    print byte_1
    aperture_width=byte_1[0]
    byte = f.read(8)
    byte_1 = numpy.frombuffer(byte, dtype=numpy.double)
    print byte_1
    aperture_height=byte_1[0]
    byte = f.read(8)
    byte_1 = numpy.frombuffer(byte, dtype=numpy.double)
    print byte_1
    delta_u=byte_1[0]
    byte = f.read(8)
    byte_1 = numpy.frombuffer(byte, dtype=numpy.double)
    print byte_1
    delta_v=byte_1[0]
    current_point = 0

    conversion_factor = (float(PSF_resolution)/float(PSF_width))
    psfs = np.zeros((20,PSF_resolution))
    while byte != "" and current_point < (sqrt_of_number_of_points*sqrt_of_number_of_points):
        x1 = numpy.frombuffer(f.read(8), dtype=numpy.double)[0] #Point three, curved surface of the lens
        y1 = numpy.frombuffer(f.read(8), dtype=numpy.double)[0]
        x2 = numpy.frombuffer(f.read(8), dtype=numpy.double)[0] #Point four, 
        y2 = numpy.frombuffer(f.read(8), dtype=numpy.double)[0]
        x3 = numpy.frombuffer(f.read(8), dtype=numpy.double)[0] #Point three, curved surface of the lens
        y3 = numpy.frombuffer(f.read(8), dtype=numpy.double)[0]
        x4 = numpy.frombuffer(f.read(8), dtype=numpy.double)[0] #Point four, 
        y4 = numpy.frombuffer(f.read(8), dtype=numpy.double)[0]
        energy_of_ray = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0]
        final_ray =  numpy.array([(x4-x3),(y4-y3)])
        refracted_vector = final_ray/np.linalg.norm(final_ray)
        intersection_point_on_curved_surface = numpy.array([x3,y3])
        scene_h_index = 0
        while scene_h_index<20:
            scene_point = numpy.array([0,(100*(scene_h_index+1))])
            scene_normal = numpy.array([0,-1])
            #scene_contact = (numpy.dot((distance_of_scene-intersection_point_on_curved_surface),(numpy.array([0,-1]))))/(numpy.dot(final_ray_normalized,numpy.array([0,-1])))
            d = numpy.dot((scene_point-intersection_point_on_curved_surface),scene_normal)/(numpy.dot(refracted_vector,scene_normal))
            scene_contact = d*refracted_vector+intersection_point_on_curved_surface
            scene_contact_bin = round(scene_contact[0]*conversion_factor + float(PSF_resolution/2))
            psfs[scene_h_index,scene_contact_bin] = psfs[scene_h_index,scene_contact_bin] + energy_of_ray
            scene_h_index = scene_h_index + 1
        current_point = current_point + 1
        #print "\draw (%f,%f) -- (%f,%f) -- (%f,%f) -- (%f,%f) [line width=0.1mm];" % (x1,y1,x2,y2,x3,y3,x4,y4)
        if current_point%10000 == 0:
            print round(float(current_point)/float(sqrt_of_number_of_points*sqrt_of_number_of_points),2)
#    recorded_array = np.zeros((height,width))
#    current_column = 0
#    current_row = 0
#    while byte != "":
#        # Do stuff with byte.
#        byte = f.read(8)
#        byte_1 = numpy.frombuffer(byte, dtype=numpy.float64)
#        if byte != "":
#            recorded_array[current_row-1,current_column-1]=byte_1
#        if current_column != width-1:
#            current_column=current_column+1
#        else:
#            current_column = 0
#            current_row = current_row+1
#        if byte_1 != 0.0:
#            print byte_1
#        if current_row == height:
#            byte = ""
#    result = recorded_array[149,:]
#    byte = f.read(4)
#    sqrt_of_number_of_points = numpy.frombuffer(byte, dtype=numpy.int32)[0]
#    print sqrt_of_number_of_points
#    byte = f.read(8)
#    lens_radius = numpy.frombuffer(byte, dtype=numpy.float64)[0]
#    print lens_radius
#    byte = f.read(8)
#    circle_height = numpy.frombuffer(byte, dtype=numpy.float64)[0]
#    print circle_height
#    current_point = 0
#    while byte != "" and current_point < (sqrt_of_number_of_points*sqrt_of_number_of_points):
#        x1 = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0] #Point 1 at the sensor
#        y1 = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0]
#        x2 = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0] #Point two, base of lens
#        y2 = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0]
#        x3 = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0] #Point three, curved surface of the lens
#        y3 = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0]
#        x4 = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0] #Point four, 
#        y4 = numpy.frombuffer(f.read(8), dtype=numpy.float64)[0]
#        if True:
#            #print "---"
#            #print x1
#            #print y1
#            #print x2
#            #print y2
#            #print x3
#            #print y3
#            #print x4
#            #print y4
#            print "\draw (%f,%f) -- (%f,%f) -- (%f,%f) -- (%f,%f) [line width=0.1mm];" % (x1,y1,x2,y2,x3,y3,x4,y4)
#        current_point = current_point + 1
#        
finally:
    f.close()

    
