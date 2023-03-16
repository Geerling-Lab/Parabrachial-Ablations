import numpy as np
import datetime
import argparse
import re
"""
Sometimes the Elitech RC-4 units get the wrong time
When you export the environmental temperatures from them, they are all off by a constant amount
This program corrects these files
As input, it requires the Elitech time and the correct time"""

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-e", "--elitech_time", help="Time from the Elitech device")
    ap.add_argument("-c", "--correct_time", help="Correct Time")
    ap.add_argument("-i", "--input_file_path", help="Path to Elitech file to be altered")
    ap.add_argument("-o", "--output_file_path", help="Path to output altered Elitech file path")
    args = vars(ap.parse_args())
    delta = datetime.datetime.strptime(args["correct_time"], "%m/%d/%y %H:%M") - \
            datetime.datetime.strptime(args["elitech_time"], "%m/%d/%y %H:%M")
    p = re.compile("\d/\d+/\\d{4} \d+:\d+:\d+ [A,P]M")
    with open(args["input_file_path"]) as input_file:
        with open(args["output_file_path"], 'w+') as output_file:
            for line in input_file:
                words = line.split("\t")
                match = p.search(line)
                if match is None:
                    print(line, end="")
                else:
                    begin, date, end = line[:match.start()], line[match.start():match.end()], line[match.end():]
                    date = datetime.datetime.strptime(date, "%m/%d/%Y %I:%M:%S %p")
                    date = date + delta
                    date = date.strftime("%m/%d/%Y %I:%M:%S %p")
                    line = begin + date + end
                    print(line, end="")
                output_file.write(line)