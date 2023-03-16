import datetime
import matplotlib.pyplot as plt
import numpy as np
import argparse
from DisplayLoggedTemperature import average_temperatures, temperatures_from_file
import glob
import os


def read_elitech_file(file_path):
    reading_header = True
    temperatures = []
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if line[:3] == "No.":
                reading_header = False
                continue
            if reading_header:
                continue
            words = line.split("      ")
            if len(words) != 3:
                continue
            _, dt, temperature = words
            try:
                dt = dt.strip()
                dt = datetime.datetime.strptime(dt, "%m/%d/%Y %I:%M:%S %p")
            except ValueError:
                print("Value Error: %s" % dt)
                continue
            temperature = float(temperature.strip())
            temperatures.append((dt, temperature))
    return temperatures


def plot_temperatures(temperatures, title):
    relative_temperatures = []
    for dt, t in temperatures:
        relative_temperatures.append(((dt - start).total_seconds() / 3600, t))
    if title == "Cage":
        xmin, xmax = ax.get_xlim()
        relative_temperatures.insert(0, (xmin, relative_temperatures[0][1]))
    relative_temperatures = np.array(relative_temperatures)
    if title == "Cage":
        line, = ax2.plot(relative_temperatures[:, 0], relative_temperatures[:, 1], "g-", label=title)
    else:
        line, = ax.plot(relative_temperatures[:, 0], relative_temperatures[:, 1], "-", label=title)
    return line


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("-d", "--directory", help="Path to directory holding temperature files")
    ap.add_argument("-t", "--datetime", help="Datetime mouse put in cold chamber")
    args = vars(ap.parse_args())

    fig, ax = plt.subplots()
    ax2 = ax.twinx()

    ax.set_ylabel("Mouse Temperature")
    ax.set_xlabel("Hours from start")
    ax.set_ylim([20, 45])
    ax2.set_ylabel("Cage Temperature")
    ax2.set_ylim([0, 45])
    start = datetime.datetime.strptime(args["datetime"], "%m/%d/%Y %I:%M:%S %p")
    lines = []
    for file in glob.glob(os.path.join(args["directory"], "*.csv")):
        title = os.path.basename(file)[:-4]
        temperatures = average_temperatures(temperatures_from_file(file))
        line = plot_temperatures(temperatures, title)
        lines.append(line)
    for file in glob.glob(os.path.join(args["directory"], "*.txt")):
        cage_temperatures = read_elitech_file(file)
        line = plot_temperatures(cage_temperatures, "Cage")
        lines.append(line)
    ax.legend(lines, [l.get_label() for l in lines])
    plt.axvline(x=0, linestyle="--", color='k')
    plt.show()