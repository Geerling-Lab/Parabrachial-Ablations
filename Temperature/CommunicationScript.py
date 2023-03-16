import serial
"""
Helper script for setting up ER-4000s
"""

ser = serial.Serial(port="COM5", timeout=1, baudrate=19200)
while True:
    s = input(">")
    ser.write(str.encode(s + "\r\n"))
    lines = ser.readlines()
    for line in lines:
        for x in line.decode("utf-8").split("\r"):
            print(x.strip())