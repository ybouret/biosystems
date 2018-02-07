import serial;

ser = serial.Serial('COM4',115200);

t_total   = 200;

line       = ser.readline();
stringline = str(line);
fields        = stringline.split();
to_write = "";
for field in fields:
    head=field.split("=")[0];
    to_write += head;
    to_write += " ";
to_write += "\n";
for field in fields:    
    value=field.split("=")[1];
    to_write += value;
    to_write += " ";

print(fields);
alpha     = fields[1].split("=")[1];
filename  = "Filtre=20_Delay=1s_alpha=%s.txt"%(alpha);
fichier   = open(filename, "w");

i = 0;

while float(fields[0].split("=")[1])<t_total:
    line       = ser.readline();
    stringline = str(line);
    fields        = stringline.split();
    
    to_write += "\n";

    for field in fields:    
        value=field.split("=")[1];
        to_write += value;
        to_write += " ";
 
    if i%10==0:
        print(fields);
                
fichier.write(to_write);

fichier.close();   
ser.close();
