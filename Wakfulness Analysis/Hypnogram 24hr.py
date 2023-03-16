import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
import datetime
from matplotlib import font_manager
import matplotlib.patches as patches
from copy import copy
import scipy.ndimage
import matplotlib.cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

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
    start_time = None
    end_time = None
    folder_location = r"R:\Fillan\Parabrachial Ablations\Sleep State\4741 Hypnogram"
    key_file = os.path.join(folder_location, "Keyfile.txt")
    with open(key_file) as f:
        start_time = f.readline().strip()
        end_time = f.readline().strip()
    start_time = datetime.datetime.strptime(start_time, "%m/%d/%Y %H:%M:%S")
    end_time = datetime.datetime.strptime(end_time, "%m/%d/%Y %H:%M:%S")

    score_file = os.path.join(folder_location, "Scores.tsv")
    band_file = os.path.join(folder_location, "Bands.tsv")
    spectra_file = os.path.join(folder_location, "Spectra.tsv")

    fig, ax = plt.subplots(nrows=4, sharex=True, gridspec_kw={'height_ratios': [1, 2, 1, 0.15]})

    scores = []
    last_rect = None
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
            hours = (dt - start_time).total_seconds() / (60 * 60)
            scores.append((hours, score))
            scores.append((hours + 0.002778, score))
            color = None
            if score == 1:
                color = (236 / 255, 250 / 255, 37 / 255, 0.2)
            elif score == 2:
                color = (0, 0, 1, 0.2)
            elif score == 3:
                color = (1, 0, 0, 0.2)
            if last_rect:
                if last_rect.get_facecolor() == color:
                    width = last_rect.get_width() + 0.002778
                    last_rect.set_width(width)
                    continue
            rect = patches.Rectangle((hours, 0), 0.002778, 3, facecolor=color, edgecolor="none")
            rects.append(rect)
            last_rect = rect
    assert scores is not None
    rect = patches.Rectangle
    scores = np.array(scores)
    ax[0].plot(scores[:, 0], scores[:, 1], 'k-', linewidth=0.5)
    ax[0].set_yticks([1, 2, 3])
    ax[0].set_yticklabels(["Wake", "NREM", "REM"])
    ax[0].tick_params(axis='both', which='major', labelsize=12)
    ax[0].tick_params(axis='both', which='minor', labelsize=8)

    rects = []
    lights_on_time = datetime.time(5, 00, 0)
    lights_off_time = datetime.time(17, 00, 0)
    current_time = start_time
    while current_time < end_time:
        start_point = (current_time - start_time).total_seconds() / 3600
        length = (end_time - current_time).total_seconds() / 3600
        if lights_on_time <= current_time.time() < lights_off_time:
            color = (255 / 255, 251 / 255, 36 / 255)
            current_time = datetime.datetime.combine(current_time.date(), lights_off_time)
        else:
            color = (30 / 255, 11 / 255, 179 / 255)
            current_time = datetime.datetime.combine(current_time.date() + datetime.timedelta(days=1), lights_on_time)
        rect = patches.Rectangle((start_point, 3.5), length, 0.5, facecolor=color)
        rects.append(rect)
    for rect in rects:
        ax[0].add_patch(copy(rect))
        ax[3].add_patch(copy(rect))
    ax[3].set_ylim([3.5, 4])
    ax[3].set_yticks([])
    ax[0].set_ylim([0.7, 4])

    times = []
    powers = []
    skip = True
    with open(spectra_file) as f:
        for line in f:
            if skip:
                if line[:4] == "Date":
                    skip = False
                continue
            date, time, _, _, *power = line.strip().split("\t")
            dt = datetime.datetime.strptime(date + " " + time, "%m/%d/%Y %H:%M:%S")
            if dt < start_time or dt > end_time:
                continue
            powers.append(list(map(float, power)))
            times.append((dt - start_time).total_seconds() / (60 * 60))
    powers = np.flip(np.array(powers).T, axis=0)
    delta = scipy.ndimage.gaussian_filter(np.sum(powers[-5:-1, :], axis=0), sigma=2)
    powers = np.log(powers)
    spectral_big = matplotlib.cm.get_cmap('nipy_spectral', 512)
    newcmp = ListedColormap(spectral_big(np.linspace(0.1, 1, 256)))
    ax[1].imshow(powers, extent=[times[0], times[-1], 0, 60], aspect='auto', cmap=newcmp, vmin=-1, vmax=6)
    ax[1].yaxis.tick_right()
    ax[1].set_ylim([0, 60])
    ax[1].set_xticks([0, 10, 20, 30, 40, 50, 60])
    ax[1].set_ylabel("EEG (Hz)")
    stem_line_width = 0.5

    """
    x = np.linspace(times[0], times[-1], num=powers.shape[1])
    _, stemlines, _ = ax[3].stem(x, delta, linefmt="k", markerfmt="none", basefmt="none")
    stemlines.set_linewidth(stem_line_width)
    ax[3].set_ylim([0, 1.2e3])
    ax[3].set_yticks([])
    ax[3].set_ylabel("Delta")
    """

    powers = []
    skip = True
    with open(band_file) as f:
        for line in f:
            if skip:
                if line[:4] == "Date":
                    skip = False
                continue
            date, time, _, _, emg, *_ = line.strip().split("\t")
            dt = datetime.datetime.strptime(date + " " + time, "%m/%d/%Y %H:%M:%S")
            if dt < start_time or dt > end_time:
                continue
            hours = (dt - start_time).total_seconds() / (60 * 60)
            emg = float(emg.strip())
            #if emg > 1800:
            #    emg = 1800
            powers.append((hours, emg))
    powers = np.array(powers)
    """
    ax[3].set_ylabel("Delta")
    delta = scipy.ndimage.gaussian_filter(powers[:, 1], sigma=1)
    _, stemlines, _ = ax[3].stem(powers[:, 0], delta, linefmt="k", markerfmt="none", basefmt="none")
    stemlines.set_linewidth(stem_line_width)
    ax[3].set_yticks([])
    """
    ax[2].set_ylabel("EMG")
    emg = scipy.ndimage.gaussian_filter(powers[:, 1], sigma=1)
    _, stemlines, _ = ax[2].stem(powers[:, 0], emg, linefmt="k", markerfmt="none", basefmt="none")
    stemlines.set_linewidth(stem_line_width)
    ax[2].set_yticks([])
    y_min, y_max = ax[2].get_ylim()
    ax[2].set_ylim(0, y_max)

    plt.xlim(powers[0, 0], powers[-1, 0])
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.xticks([0, 4, 8, 12, 16, 20, 24])
    plt.xlabel("Zeitgeber Time")
    plt.show()

