# -*- coding: utf-8 -*-
#______________________________________________________________________________
#
# used modules
#______________________________________________________________________________
import serial;
import time;

#______________________________________________________________________________
#
# create arduino as a serial interface
#______________________________________________________________________________

print("---> Initialize serial com");
serial_bauds = 115200;
serial_path  = '/dev/cu.usbmodem1421';
arduino     = serial.Serial(port=serial_path,baudrate=serial_bauds);
print("---> Cleanup Arduino");
arduino.setDTR(False);
time.sleep(1);
arduino.flushInput();
arduino.setDTR(True);

#______________________________________________________________________________
#
# global variables
#______________________________________________________________________________
fields   = [];
t        = float(0);
warmup   = 5;
record   = 10;
filename = "output.dat";

#______________________________________________________________________________
#
# SendCommand function
#______________________________________________________________________________
def SendCommand(cmd):
    global arduino;
    output = (str(cmd) );
    print("---> Sending '%s'"%output);
    output = output + "\n";
    arduino.write( output.encode() );
    arduino.flush();
    return;

#______________________________________________________________________________
#
# ReadFields function
#______________________________________________________________________________
def ReadFields():
    global fields;
    global t;
    global arduino;
    fields = [];
    line = str(arduino.readline()).strip();
    print("arduino:<%s>"%line);
    if len(line) <= 0:
        return False; # empty line
    if '#' == line[0]:
        return False; # comment like output, discard
    fields = line.split();
    if len(fields) <= 2:
        return False; # at least t=xxx angle=xxx
    t_info = fields[0].split('=');
    if len(t_info) != 2:
        return False; # field is not 'field=value'
    if  't' != t_info[0]:
        return False; # fisrt field is 't=xxx'
    t = float(t_info[1]);
    return True;
    
#______________________________________________________________________________
#
# saving data
#______________________________________________________________________________   
headers  = False; 

def SaveData() :
    global headers;
    global filename;
    global fields;
    global t;
    output_range= range(1,len(fields));
    #__________________________________________________________________________
    #
    # one time write fields names
    #__________________________________________________________________________
    if(not headers):
        #write fields name, first one is time
        with open(filename,"w") as fp:
            fp.write('#t');
            for f in output_range:
                fp.write(' ');
                fp.write(fields[f].split('=')[0]);
            fp.write('\n');
        headers=True;
   
    #__________________________________________________________________________
    #
    # write fields value
    #__________________________________________________________________________
    with open(filename,"a") as fp:
        fp.write(str(t));
        for f in output_range:
            fp.write(' ');
            fp.write(fields[f].split('=')[1]);
        fp.write('\n');
    return;
    
#______________________________________________________________________________
#
# run: SEND A COMMENT FROM ARDUINO IN SETUP!!!
#______________________________________________________________________________  
print("---> Running");
ReadFields();
SendCommand("period 4");
SendCommand("motion tri");
run_time = warmup+record;
while arduino.is_open:
    if not ReadFields():
        continue;
        
    if(t<=warmup):
        continue;
        
    SaveData();
    
    if t>=run_time: 
        break;
print("---> done");

SendCommand("period 10");
ReadFields();
SendCommand("amplitude 30");
ReadFields();
SendCommand("motion tri");
ReadFields();

arduino.close();

