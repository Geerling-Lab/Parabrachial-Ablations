import numpy as np
import matplotlib.pyplot as plt
import argparse
import datetime
import collections
from matplotlib import font_manager
"""
This program takes in two files;
A file containing the number of movements per epoch and a file containing the EEG scoring
"""

if __name__ == "__main__":
    font_path = r"C:\Users\Fillan\Documents\Fonts\Myriad Pro Regular.ttf"
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = prop.get_name()
    plt.rcParams['font.size'] = 18
    plt.rcParams['axes.linewidth'] = 2

    ap = argparse.ArgumentParser()
    ap.add_argument("-m", "--movement", help="Path to movement file from VideoScoring.py")
    ap.add_argument("-e", "--eeg", help="Path to EEG file exported from Sirenia")
    args = vars(ap.parse_args())

    times = []
    scores = []
    movements = []
    start_time = None
    with open(args["eeg"]) as f:
        header = True
        for line in f:
            if header:
                if line[:4] == "Date":
                    header = False
                continue
            date, dt_time, _, _, score = line.strip().split("\t")
            epoch_dt = datetime.datetime.strptime(date + " " + dt_time.strip(), "%m/%d/%Y %H:%M:%S")
            if start_time is None:
                start_time = epoch_dt
            time_since_start = (epoch_dt - start_time).total_seconds()
            times.append(time_since_start)
            scores.append(int(float(score)))

    dt_counter = 0
    with open(args["movement"]) as f:
        for line in f:
            epoch_str, movement = line.strip().split("\t")
            movement = int(movement)
            if movement > 75000:
                movement = 75000
            epoch_dt = datetime.datetime.strptime(epoch_str, "%m/%d/%Y %H:%M:%S")
            time_since_start = (epoch_dt - start_time).total_seconds()
            while time_since_start > times[dt_counter]:
                movements.append(movement)
                dt_counter += 1
                if dt_counter >= len(times):
                    break
            if dt_counter >= len(times):
                break

    times = np.array(times)
    scores = np.array(scores)
    movements = np.array(movements)

    fig, ax = plt.subplots(2, sharex=True)
    ax[0].set_yticks([1, 2, 3])
    ax[0].set_yticklabels(["Wake", "NREM", "REM"])
    ax[0].plot(times, scores, "k-")
    ax[1].plot(times, movements, "k-")
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.show()

    fig, ax = plt.subplots()
    i = 0
    names = {1: "Wake", 2: "NREM", 3: "REM"}
    ax.boxplot([movements[scores==1], movements[scores==2], movements[scores==3]], 0, "")
    ax.set_xticklabels(["Wake", "NREM", "REM"])
    ax.set_ylabel("Movement")
    #ax.set_yscale("log")
    plt.show()

    fig, ax = plt.subplots(nrows=3, sharex=True)
    for i in [0, 1, 2]:
        logbins = np.geomspace(1e2, 1e5, 50)
        bins = np.arange(0, 2e5, 2e3)
        ax[i].hist(movements[scores==i+1], bins=logbins)
        ax[i].set_title(i)
        ax[i].set_xscale('log')

    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.xlabel("Movement")
    plt.ylabel("Number of pixels")
    plt.show()

