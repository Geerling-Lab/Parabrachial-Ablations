"""
Helper script for editing out stray E-mitter recordings.
Some of the ER-4000s have a poor signal, and will read the E-mitter on the neighboring ER-4000. For example,
ER-4000 #4 will read the E-mitter on ER-4000 #3. These can be distinguished because they have unusually
high or low temperatures, because of the different calibration.

This script reads through a .csv file output by "TemperatureLogger.py", and removes any lines that are above a
threshold maximum or below a threshold minimum temperature
"""
input_file_path = r"R:\Fillan\Parabrachial Ablations\Temperature\6314 6315 6316 6317\Control 11 1 22\Core 6317 OLD.csv"
output_file_path = r"R:\Fillan\Parabrachial Ablations\Temperature\6314 6315 6316 6317\Control 11 1 22\Core 6317.csv"
maximum_temperature = 37.5
minimum_temperature = 0
counter = 0
with open(input_file_path) as i_f:
    with open(output_file_path, "w+") as o_f:
        for line in i_f:
            dt, temperature = line.strip().split(",")
            temperature = float(temperature)
            if minimum_temperature < temperature < maximum_temperature or temperature == 0:
                o_f.write("%s,%.4f\n" % (dt, temperature))
            else:
                counter += 1
print("%i removed" % counter)

