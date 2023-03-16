import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import datetime
"""
Generates a figure showing raw EEG and EMG signals during stimulation
"""


def read_raw_file(raw_file_path):
    df = pd.read_csv(raw_file_path, sep="\t", header=5)
    start_time = datetime.datetime.strptime(df["Date"][0].strip(), "%H:%M:%S")
    end_time = datetime.datetime.strptime(df["Date"][-1].strip(), "%H:%M:%S")
    eeg = np.array(df['Time from Start'])
    emg = np.array(df['EEG2'])
    times = np.array(df['Time Stamp'])
    return start_time, end_time, eeg, emg, times


def read_annotation_file(annotation_file_path, start_time, end_time):
    df = pd.read_csv(annotation_file_path, sep="\t", header=4)
    annotation_times = list(df.loc[df["Channel"] == "EEG2"]["Start Time"])
    valid_annotation_times = []
    for time in annotation_times:
        time = datetime.datetime.strptime(time.split(" ")[1], "%H:%M:%S.%f")
        if start_time < time < end_time:
            valid_annotation_times.append((time - start_time).total_seconds())
    return valid_annotation_times


if __name__ == "__main__":
    raw_file_path = r"R:\EEG Scores\Optogenetic EEG Scores\Raw\JR24 Stim.tsv"
    annotation_file_path = "R:\EEG Scores\Optogenetic EEG Scores\Annnotations\JR24 Annotations.tsv"
    start_time, end_time, eeg, emg, times = read_raw_file(raw_file_path)
    annotations = read_annotation_file(annotation_file_path, start_time, end_time)
    first_annotation = annotations[0]
    fig, ax = plt.subplots(nrows=3, sharex=True, gridspec_kw={'height_ratios': [5, 5, 1]})
    plt.subplots_adjust(wspace=0, hspace=-0.01)
    ax[0].margins(x=0, y=0)
    ax[1].margins(x=0, y=0)
    ax[2].margins(x=0, y=0)
    ax[0].set_yticks([0])
    ax[0].set_yticklabels(["EEG"])
    ax[1].set_yticks([0])
    ax[1].set_yticklabels(["EMG"])
    ax[2].set_yticks([])
    ax[0].spines['bottom'].set_visible(False)
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['bottom'].set_visible(False)
    ax[2].spines['top'].set_visible(False)
    ax[0].plot(times - first_annotation, eeg, lw=0.3, c='k')
    ax[1].plot(times - first_annotation, emg, lw=0.3, c='k')

    for annotation in annotations:
        rect = patches.Rectangle((annotation - first_annotation, 0), 0.02, 1, linewidth=0, edgecolor='none', facecolor="#1FF4FF")
        ax[2].add_patch(rect)
    ax[2].set_ylim([-1, 1])
    ax[2].set_xlim([-30, 40])
    plt.show()
