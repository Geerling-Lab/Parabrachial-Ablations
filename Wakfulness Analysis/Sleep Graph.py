import numpy as np
import matplotlib.pyplot as plt
import datetime
import re
import argparse
import glob
import os
import random
import pandas as pd
from enum import Enum
from itertools import repeat
import Helper

import scipy.stats
from matplotlib import font_manager
import matplotlib.patches as patches


class Vector(Enum):
    Cas3 = 1
    mCherry = 2


class Mouse:
    def __init__(self, number, vector, ZT, stereological_count, movement):
        self.times = np.zeros(shape=(25, 3), dtype=np.int16)
        self.fraction_times = np.zeros(shape=self.times.shape)
        self.number = number
        self.vector = vector
        self.color = 'r'
        if self.vector == Vector.mCherry:
            self.color = "gray"
        self.ZT = ZT
        self.valid_scoring = False
        self.n = 1
        self.stereological_count = stereological_count
        self.total_wake = np.nan
        self.movement = movement
        self.last_states = []
        self.bout_length = 0
        self.bout_state = None
        self.bouts = [[], [], []]
        self.transition_matrix = np.zeros(shape=(3, 3), dtype=np.int32)

    def add_time(self, hour, state):
        if state > 3:
            return
        self.last_states.insert(0, state)
        while len(self.last_states) > 5:
            del self.last_states[5]
        self.valid_scoring = True
        self.times[hour, self.last_states[0] - 1] += 1
        if len(self.last_states) > 1:
            self.transition_matrix[self.last_states[1] - 1, self.last_states[0] - 1] += 1
            if self.bout_state is None:
                self.bout_state = self.last_states[0]
                self.bout_length = 0
            if self.last_states[0] != self.bout_state and self.last_states[1] != self.bout_state:
                if self.last_states[0] == self.last_states[1] or self.last_states[2] != self.bout_state:
                    if self.bout_length > 1:
                        self.bouts[self.bout_state - 1].append(self.bout_length)
                    self.bout_state = self.last_states[0]
                    self.bout_length = 1
            self.bout_length += 1

    def finish_times(self):
        self.times[24, :] = self.times[0, :]
        if np.sum(self.times[:24, :]) > 10000: # If longer than 10,000 epochs of data was collected, then it's two days
            self.n = 2
            self.transition_matrix = (self.transition_matrix / 2).astype(np.int32)
        for i in range(3):
            self.fraction_times[:, i] = self.times[:, i] / np.sum(self.times, axis=1)
        self.fraction_times[self.fraction_times == np.nan] = 0
        self.total_wake = np.nansum(self.times[:, 0]) / np.nansum(self.times)

    def plot(self, axs):
        if self.valid_scoring:
            for i in range(3):
                axs[i].plot(np.arange(25), self.fraction_times[:, i], color=self.color)

    def plot_scatter(self, ax):
        print(self.stereological_count, self.total_wake, self.color)
        ax.plot(self.stereological_count, self.total_wake, color=self.color)

    def get_bouts(self, i, j):
        # j is sort of a dummy variable
        # When j = 0, return the number of bouts
        # When j = 1, return the average length of bouts
        if j == 0:
            return len(self.bouts[i])
        elif j == 1:
            return sum(self.bouts[i]) / len(self.bouts[i])
        else:
            raise AssertionError("Invalid argument j=%s" % j)


def error_bar(ax, x, y, yerr, color, label=None):
    ax.plot(x, y, color, label=label)
    lower = [i - j for i, j in zip(y, yerr)]
    upper = [i + j for i, j in zip(y, yerr)]
    ax.fill_between(x, lower, upper, color=color, alpha=0.2, edgecolor='w')


def graph_versus_stereology(mice, func, ylabel):
    fig, ax = plt.subplots()
    Cas3_stereological_counts, Cas3_y_values = [], []
    mCh_stereological_counts, mCh_y_values = [], []
    for mouse_number, mouse in mice.items():
        y_value = func(mouse)
        if not np.isnan(y_value) and not np.isnan(mouse.stereological_count):
            if mouse.vector == Vector.mCherry:
                mCh_stereological_counts.append(mouse.stereological_count)
                mCh_y_values.append(y_value)
            else:
                Cas3_stereological_counts.append(mouse.stereological_count)
                Cas3_y_values.append(y_value)
    ax.scatter(mCh_stereological_counts, mCh_y_values, c="gray")
    ax.scatter(Cas3_stereological_counts, Cas3_y_values, c="r")
    _, p_pearson = scipy.stats.pearsonr(mCh_stereological_counts + Cas3_stereological_counts, mCh_y_values + Cas3_y_values)
    print("********")
    print(ylabel)
    print("********")
    print("Correlation: p=%.3f" % p_pearson)
    mCh_y_values = np.array(mCh_y_values)
    Cas3_y_values = np.array(Cas3_y_values)
    print("mCh mean: %.2f ± %.2f (n=%i)" % (np.mean(mCh_y_values), np.std(mCh_y_values), len(mCh_y_values)))
    print("Cas3 mean: %.2f ± %.2f (n=%i)" % (np.mean(Cas3_y_values), np.std(Cas3_y_values), len(Cas3_y_values)))
    _, p_ttest = scipy.stats.ttest_ind(Cas3_y_values, mCh_y_values)
    print("T-test: p=%.3f" % p_ttest)
    plt.plot([0, 0], [.5, .5], 'w-', label="p=%.3f" % p_pearson)
    plt.legend()
    plt.xlabel("L10GFP (Vglut2) neurons remaining")
    plt.ylabel(ylabel)
    plt.ylim([0, 1])
    plt.show()


if __name__ == "__main__":
    font_path = r"C:\Users\Fillan\Documents\Fonts\Myriad Pro Regular.ttf"
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = prop.get_name()
    plt.rcParams['font.size'] = 18
    plt.rcParams['axes.linewidth'] = 2

    parser = argparse.ArgumentParser(description="Creates a Hypnogram from a Sirennia TSV file")
    parser.add_argument("-d", "--directory", type=str, help="Directory with tsv files")
    parser.add_argument("-o", "--output", type=str, help="Output filepath")
    args = parser.parse_args()
    mice = {}
    display_individual = True
    display_average = False
    dfs = pd.read_excel(io=r"R:\Fillan\Parabrachial Ablations\Lesion Mice Overview.xlsx", sheet_name="Overview")
    for index, row in dfs.iterrows():
        if "mCh" in row["Vector"]:
            vector = Vector.mCherry
        elif "Cas3" in row["Vector"]:
            vector = Vector.Cas3
        else:
            continue
        if not isinstance(row["Light on"], datetime.time):
            pass
        ZT = row["Light on"]
        stereological_count = row["Stereology"]
        movement = row["Movement"] / 100
        number = "%05d" % row["Mouse"]
        mice[number] = Mouse(number=number, vector=vector, ZT=ZT, stereological_count=stereological_count, movement=movement)
    for name in glob.glob('%s/*.tsv' % args.directory):
        digits = re.findall("\d{5}", name)
        if len(digits) > 0:
            number = digits[0]
            print(number)
            if number not in mice:
                raise IOError("%s not in Lesion Overview.xlsx file" % number)
            with open(name) as f:
                for line in f:
                    words = line.strip().split("\t")
                    if len(words) != 5:
                        continue
                    if words[0].count("/") != 2:
                        continue
                    datewithtime = datetime.datetime.strptime(words[0] + "," + words[1].strip(), "%m/%d/%Y,%H:%M:%S")
                    date = datewithtime.date()
                    time = datewithtime.time()
                    time_from_ZT = datewithtime - datetime.datetime.combine(date, mice[number].ZT)
                    if time_from_ZT < datetime.timedelta(0):
                        time_from_ZT += datetime.timedelta(1, 0, 0)
                    hour = time_from_ZT.seconds // 3600
                    mice[number].add_time(hour, int(float(words[4])))
            mice[number].finish_times()
    fig, axs = plt.subplots(nrows=3, ncols=3)
    labels = ["Wake", "NREM", "REM"]
    for i in range(3):
        for j in range(3):
            Cas3_transitions = []
            mCh_transitions = []
            for mouse_number, mouse in mice.items():
                if mouse.valid_scoring:
                    if mouse.vector == Vector.Cas3:
                        Cas3_transitions.append(mouse.transition_matrix[i, j])
                    else:
                        mCh_transitions.append(mouse.transition_matrix[i, j])
            _, p = scipy.stats.ttest_ind(mCh_transitions, Cas3_transitions)
            axs[i, j].set_title("%s -> %s, p=%.3f" % (labels[i], labels[j], p))
            axs[i, j].scatter(Helper.fillan_bee_swarm(mCh_transitions, 1, 0.5, 0.05), mCh_transitions, color="lightgrey")
            mCh_average = sum(mCh_transitions) / len(mCh_transitions)
            axs[i, j].plot([0.8, 1.2], [mCh_average, mCh_average], color="k")
            axs[i, j].scatter(Helper.fillan_bee_swarm(Cas3_transitions, 2, 0.5, 0.05), Cas3_transitions, color="r")
            Cas3_average = sum(Cas3_transitions) / len(Cas3_transitions)
            axs[i, j].plot([1.8, 2.2], [Cas3_average, Cas3_average], color="k")
            axs[i, j].set_xlim([0.5, 2.5])
            axs[i, j].set_xticks([1, 2])
            _, y_max = axs[i, j].get_ylim()
            axs[i, j].set_ylim([0, y_max])
            axs[i, j].set_xticklabels(["mCh", "Cas3"])
            print("%s -> %s, p=%.3f" % (labels[i], labels[j], p))
            print("mCh: %.2f ± %.2f" % (mCh_average, np.std(mCh_transitions)))
            print("Cas3: %.2f ± %.2f" % (Cas3_average, np.std(Cas3_transitions)))
    plt.show()
    fig, axs = plt.subplots(nrows=2, ncols=3)
    functions = ["Number", "Average Length"]
    for i in range(3):
        for j in range(2):
            Cas3_bouts, mCh_bouts = [], []
            for mouse_number, mouse in mice.items():
                if mouse.valid_scoring:
                    if mouse.vector == Vector.Cas3:
                        Cas3_bouts.append(mouse.get_bouts(i, j))
                    else:
                        mCh_bouts.append(mouse.get_bouts(i, j))
            print("******")
            print("%s %s bouts" % (functions[j], labels[i]))
            print("******")
            mCh_bouts = np.array(mCh_bouts)
            Cas3_bouts = np.array(Cas3_bouts)
            print("mCh mean: %.3f ± %.3f (n=%i)" % (np.mean(mCh_bouts), np.std(mCh_bouts), np.size(mCh_bouts)))
            print("Cas3 mean: %.3f ± %.3f (n=%i)" % (np.mean(Cas3_bouts), np.std(Cas3_bouts), np.size(Cas3_bouts)))
            axs[j, i].plot([0.8, 1.2], [np.mean(mCh_bouts), np.mean(mCh_bouts)], "k-")
            axs[j, i].plot([1.8, 2.2], [np.mean(mCh_bouts), np.mean(mCh_bouts)], "k-")
            axs[j, i].scatter(Helper.fillan_bee_swarm(mCh_bouts, 1, 0.5, 0.05), mCh_bouts, color="lightgrey")
            axs[j, i].scatter(Helper.fillan_bee_swarm(Cas3_bouts, 2, 0.5, 0.05), Cas3_bouts, color="r")
            _, p = scipy.stats.ttest_ind(Cas3_bouts, mCh_bouts)
            axs[j, i].set_title("%s %s bouts: p=%.3f" % (functions[j], labels[i], p))
            axs[j, i].set_xlim([0.5, 2.5])
            axs[j, i].set_xticks([1, 2])
            axs[j, i].set_xticklabels(["mCh", "Cas3"])
            if j == 1:
                axs[j, i].set_ylabel("Epochs")
            _, y_max = axs[j, i].get_ylim()
            axs[j, i].set_ylim([0, y_max])
    plt.tight_layout()
    plt.show()
    Cas3 = []
    Cas3_n = 0
    mCh = []
    mCh_n = 0
    fig, axs = plt.subplots(nrows=3, sharex=True)
    plt.subplots_adjust(wspace=0, hspace=0)
    for mouse_number, mouse in mice.items():
        if mouse.valid_scoring:
            mouse.plot(axs)
            if mouse.vector is Vector.Cas3:
                Cas3_n += 1
                Cas3.append(mouse.fraction_times)
                if mouse.n == 2:
                    Cas3.append(mouse.fraction_times)
            if mouse.vector is Vector.mCherry:
                mCh_n += 1
                mCh.append(mouse.fraction_times)
                if mouse.n == 2:
                    mCh.append(mouse.fraction_times)
    for i in range(3):
        axs[i].set_ylim([0, 1])
    axs[0].set_ylabel("Wake")
    axs[1].set_ylabel("NREM")
    axs[2].set_ylabel("REM")
    plt.xticks([0, 4, 8, 12, 16, 20, 24])
    plt.xlabel("Zeitgeber Time")
    plt.xlim([0, 24])
    plt.show()

    Cas3 = np.array(Cas3)
    mCh = np.array(mCh)
    Cas3_mean = np.mean(Cas3, axis=0)
    mCh_mean = np.mean(mCh, axis=0)
    Cas3_stddev = np.std(Cas3, axis=0)
    mCh_stddev = np.std(mCh, axis=0)
    fig, axs = plt.subplots(nrows=4, sharex=True, gridspec_kw={'height_ratios': [1, 1, 1, 0.1]})
    rect = patches.Rectangle((0, 0), 12, 1, facecolor=(255 / 255, 251 / 255, 36 / 255))
    axs[3].add_patch(rect)
    rect = patches.Rectangle((12, 0), 12, 1, facecolor=(30 / 255, 11 / 255, 179 / 255))
    axs[3].add_patch(rect)
    axs[3].set_ylim([0, 1])
    axs[3].set_yticks([])
    Cas3_total = np.sum(Cas3[:, :24, :], axis=1)
    mCh_total = np.sum(mCh[:, :24, :], axis=1)
    plt.subplots_adjust(wspace=0, hspace=0)
    for i in range(3):
        print("****")
        print(i)
        print("****")
        print("Cas3 mean: %.3f ± %.3f" % (np.mean(Cas3_total[:, i]), np.std(Cas3_total[:, i])))
        print("mCh mean: %.3f ± %.3f" % (np.mean(mCh_total[:, i]), np.std(mCh_total[:, i])))
        _, p = scipy.stats.ttest_ind(Cas3_total[:, i], mCh_total[:, i])
        print("Ttest: p=%.3f" % p)
        error_bar(ax=axs[i], x=np.arange(25), y=mCh_mean[:, i], yerr=mCh_stddev[:, i], color='gray', label="mCh (n=%i)" % mCh_n)
        error_bar(ax=axs[i], x=np.arange(25), y=Cas3_mean[:, i], yerr=Cas3_stddev[:, i], color='r', label="Cas3 (n=%i)" % Cas3_n)
    for i in range(3):
        axs[i].set_ylim([0, 1])
        axs[i].set_yticks([0, 1])
        axs[i].set_xticks([])
        axs[i].tick_params(axis='y', which='major', labelsize=10)
    axs[0].legend(prop={'size': 8})
    axs[0].set_ylabel("Wake")
    axs[1].set_ylabel("NREM")
    axs[2].set_ylabel("REM")
    axs[3].set_xticks([0, 4, 8, 12, 16, 20, 24])
    plt.xlabel("Zeitgeber Time")
    plt.xlim([0, 24])
    plt.show()

    graph_versus_stereology(mice, func=lambda x: x.movement, ylabel="Activity")
    graph_versus_stereology(mice, func=lambda x: x.total_wake, ylabel="Total wakefulness")
    graph_versus_stereology(mice, func=lambda x: np.nansum(x.times[:12, 0]) / np.nansum(x.times[:12, :]), ylabel="Light phase wakefulness")
    graph_versus_stereology(mice, func=lambda x: np.nansum(x.times[12:24, 0]) / np.nansum(x.times[12:24, :]), ylabel="Dark phase wakefulness")