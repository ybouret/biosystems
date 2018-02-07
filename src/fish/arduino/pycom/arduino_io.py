# -*- coding: utf-8 -*-

import serial;
import time;

print("---> initialize serial com");
serial_bauds = 115200;
serial_path  = '/dev/cu.usbmodem21';
serialIO     = serial.Serial(port=serial_path,baudrate=serial_bauds);
print("---> waiting");
time.sleep(0.5);

print("---> sending parameters");
serialIO.write("period 5\n".encode());
serialIO.flush();

print("---> reading data");

#______________________________________________________________________________
#
#get a valid line and a start time
#______________________________________________________________________________
time_start = 0;
while serialIO.is_open:
    # read a string line, stripping trailing whitespaces, and split fields
    fields = str(serialIO.readline()).strip().split();
    if len(fields)<=2:
        continue;
    #the first field is time info
    time_field = fields[0];
    time_info  = time_field.split('=');
    if len(time_info) != 2:
        continue; # bad format...
    if time_info[0] != 't':
        continue; # bad format...
    time_start = float(time_info[1]);
    break;
print("time_start=%g"%time_start);

#______________________________________________________________________________
#
#read other lines, assuming well formatted
#______________________________________________________________________________
while True:
     # read a string line, stripping trailing whitespaces
    fields     = str(serialIO.readline()).strip().split();
    num_fields = len(fields);
    if num_fields<=2:
        print("bad number of fields");
        break; #error
    time_info = fields[0].split('=');
    if len(time_info) != 2:
        print("bad time format");
        break; #error
    if time_info[0] != 't':
        print("first field is not 't'");
        break; #error
    t = float(time_info[1]);
    if t-time_start>5:
        break; # normal ending
    # process fields for saving
    print(fields);
    
print("---> done");
serialIO.close();
