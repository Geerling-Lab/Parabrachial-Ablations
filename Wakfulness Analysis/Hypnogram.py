import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
import datetime
from matplotlib import font_manager
import matplotlib.patches as patches
from copy import copy

"""
Generates a pretty hypnogram for inclusion in a paper
Requires two files:
  Score file, exported from Sirenia
  Raw file, with EEG2, EMG, EEG2 delta power, and EEG2 theta power
This also requires an start and end time for the hypnogram, set with flags
    -s: start time, in format "dd-mm-yyyy hh:mm:ss
    -e: end time, in format "dd-mm-yyyy hh:mm:ss
This will produce a matplotlib graph with the 
scores, raw waveforms of EEG2, EMG, and theta and delta power displayed"""

if __name__ == "__main__":
    font_path = r"C:\Users\Fillan\Documents\Fonts\Myriad Pro Regular.ttf"
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = prop.get_name()
    plt.rcParams['font.size'] = 18
    plt.rcParams['axes.linewidth'] = 2

    ap = argparse.ArgumentParser()
    ap.add_argument("-s", "--start", help="Hypnogram start time, in format dd-mm-yyyy hh:mm:ss")
    ap.add_argument("-e", "--end", help="Hypnogram end time, in format dd-mm-yyyy hh:mm:ss")
    args = vars(ap.parse_args())
    start_time = datetime.datetime.strptime(args["start"], "%m-%d-%Y %H:%M:%S")
    end_time = datetime.datetime.strptime(args["end"], "%m-%d-%Y %H:%M:%S")
    folder_location = r"C:\Users\Filla\Downloads\4739 Hypnogram"
    score_file = os.path.join(folder_location, "Scores.tsv")
    raw_file = os.path.join(folder_location, "Raw.txt")

    scores = []
    fig, ax = plt.subplots(nrows=5, sharex=True)
    rects = []
    with open(score_file) as f:
        skip = True
        for line in f:
            if skip:
                if line[:4] == "Date":
                    skip = False
                continue
            date, time, _, _, score = line.strip().split("\t")
            dt = datetime.datetime.strptime(date + " " + time, "%m/%d/%Y %H:%M:%S")
            if dt < start_time or dt > end_time:
                continue
            score = int(float(score.strip()) + 0.5)
            minutes = (dt - start_time).total_seconds() / 60
            scores.append((minutes, score))
            scores.append((minutes + 0.16666, score))
            color = None
            if score == 1:
                color = (236 / 255, 250 / 255, 37 / 255, 0.2)
            elif score == 2:
                color = (0, 0, 1, 0.2)
            elif score == 3:
                color = (1, 0, 0, 0.2)
            rects.append(patches.Rectangle((minutes, -1e3), 0.16667, 2e4, facecolor=color, edgecolor="none"))
    scores = np.array(scores)
    ax[4].plot(scores[:, 0], scores[:, 1], 'k-')
    ax[4].set_yticks([1, 2, 3])
    ax[4].set_yticklabels(["Wake", "NREM", "REM"])
    ax[4].tick_params(axis='both', which='major', labelsize=12)
    ax[4].tick_params(axis='both', which='minor', labelsize=8)

    raws = []
    powers = [[], [], []]
    record_start = None
    with open(raw_file) as f:
        skip = True
        for line in f:
            if skip:
                if line[:4] == "Date":
                    skip = False
                continue
            score_index = 0
            date, time, _, seconds_from_start, eeg, emg, *extra = line.strip().split("\t")
            dt = datetime.datetime.strptime(date + " " + time, "%m/%d/%Y %H:%M:%S")
            if not record_start:
                record_start = (start_time - dt).total_seconds()
            if dt < start_time or dt > end_time:
                continue
            eeg = float(eeg)
            emg = float(emg)
            minutes = (float(seconds_from_start) - record_start) / 60
            raws.append((minutes, eeg, emg))
            if extra:
                theta, delta = extra
                delta = int(float(delta) + 0.5)
                theta = int(float(theta) + 0.5)
                while scores[score_index + 1][0] <= minutes:
                    score_index += 1
                powers[int(scores[score_index][1] - 0.5)].append((minutes + 0.0833, delta, theta))
                """
                Add 0.8333 (5 seconds in minutes, or one half of an epoch) to make things line up
                """
    raws = np.array(raws)
    ax[0].plot(raws[:, 0], raws[:, 1], "k-", linewidth=0.5)
    ax[0].set_yticks([])
    ax[0].set_ylabel("EEG")
    ax[1].plot(raws[:, 0], raws[:, 2], "k-", linewidth=0.5)
    ax[1].set_ylabel("EMG")
    ax[1].set_yticks([])
    wake_powers = np.array(powers[0])
    nrem_powers = np.array(powers[1])
    rem_powers = np.array(powers[2])
    stem_line_width = 3
    if wake_powers.size > 0:
        _, stemlines, _ = ax[2].stem(wake_powers[:, 0], wake_powers[:, 1], linefmt="k", markerfmt="none", basefmt="none")
        stemlines.set_linewidth(stem_line_width)
    if nrem_powers.size > 0:
        _, stemlines, _ = ax[2].stem(nrem_powers[:, 0], nrem_powers[:, 1], linefmt="k", markerfmt="none", basefmt="none")
        stemlines.set_linewidth(stem_line_width)
    if rem_powers.size > 0:
        _, stemlines, _ = ax[2].stem(rem_powers[:, 0], rem_powers[:, 1], linefmt="k", markerfmt="none", basefmt="none")
        stemlines.set_linewidth(stem_line_width)
    ax[2].set_ylabel("Delta")
    ax[2].set_yticks([])
    if wake_powers.size > 0:
        _, stemlines, _ = ax[3].stem(wake_powers[:, 0], wake_powers[:, 2], linefmt="k", markerfmt="none", basefmt="none")
        stemlines.set_linewidth(stem_line_width)
    if nrem_powers.size > 0:
        _, stemlines, _ = ax[3].stem(nrem_powers[:, 0], nrem_powers[:, 2], linefmt="k", markerfmt="none", basefmt="none")
        stemlines.set_linewidth(stem_line_width)
    if rem_powers.size > 0:
        _, stemlines, _ = ax[3].stem(rem_powers[:, 0], rem_powers[:, 2], linefmt="k", markerfmt="none", basefmt="none")
        stemlines.set_linewidth(stem_line_width)
    ax[3].set_ylabel("Theta")
    ax[3].set_yticks([])

    for i in range(5):
        y_min, y_max = ax[i].get_ylim()
        for rect in rects:
            ax[i].add_patch(copy(rect))
        if i == 4:
            y_max = 3.5
        ax[i].set_ylim([y_min, y_max])
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xlabel("Minutes")
    plt.show()

