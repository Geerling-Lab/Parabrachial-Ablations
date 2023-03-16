import numpy as np
import matplotlib.pyplot as plt
import glob
import os
from matplotlib.patches import Rectangle
from matplotlib import font_manager
import datetime
import argparse
import scipy.stats

"""
This program analyzes the data from the "heat lamp" experiment. I connected a nose poke to a Raspberry Pi. Every time the mouse pokes his nose in the port, a heat lamp comes on for ten seconds,
and the nose poke is logged. The Raspberry Pi logs the data in the following format
85-86.5
100.9-115.1
132.1-133.3
Where each line is the beginning and end of a nosepoke, as measured in seconds from the start of the experiment.
The first nose poke is considered the start of the experiment.

These files should be saved together in a folder, with the mouse's number first, followed by the date. For example,
6000 2 27 22 10 00.txt <- Mouse 6000, trial started at 10:00 AM on 2/27/22
6000 2 28 22 10 00.txt
6000 2 29 22 10 00.txt
6001 2 27 22 12 00.txt
6001 2 28 22 12 00.txt
"""
TOTAL_DAYS = 7


class OneDayData:
    def __init__(self):
        self.null = False
        self.pokes = [0, 0]
        self.time_on = [0, 0]

    def __iadd__(self, other):
        assert type(other) == OneDayData
        self.pokes[0] += other.pokes[0]
        self.pokes[1] += other.pokes[1]
        self.time_on[0] += other.time_on[0]
        self.time_on[1] += other.time_on[1]
        return self

    def __truediv__(self, other):
        self.pokes[0] /= other
        self.pokes[1] /= other
        self.time_on[0] /= other
        self.time_on[1] /= other
        return self


class MouseData:
    def __init__(self, mouse, vector):
        self.mouse = mouse
        self.vector = vector
        self.values = []
        self.average_cold_pokes = None
        self.average_hot_pokes = None
        for _ in range(TOTAL_DAYS):
            self.values.append(OneDayData())

    def input(self, day_number, channel, total_pokes, total_time_on):
        if total_pokes is None or total_time_on is None:
            self.values[day_number].null = True
            return
        if channel == "1" or channel == "3":
            self.values[day_number].pokes[0] += total_pokes
            self.values[day_number].time_on[0] += total_time_on
        else:
            self.values[day_number].pokes[1] += total_pokes
            self.values[day_number].time_on[1] += total_time_on

    def calculate_average_pokes(self):
        self.average_cold_pokes = OneDayData()
        number_cold_days = 0
        for value in [self.values[4]]:#self.values[1:-1]:
            if not value.null:
                self.average_cold_pokes += value
                number_cold_days += 1
        if number_cold_days >= 1:
            self.average_cold_pokes = self.average_cold_pokes / number_cold_days

        self.average_hot_pokes = OneDayData()
        number_hot_days = 0
        for value in [self.values[0]]:#, self.values[-1]]:
            if not value.null:
                self.average_hot_pokes += value
                number_hot_days += 1
        if number_hot_days:
            self.average_hot_pokes = self.average_hot_pokes / number_hot_days


def file_name_key(a):
    a = os.path.basename(a)
    a = " ".join(a[:-6].split(" ")[1:])
    datetime_a = datetime.datetime.strptime(a, "%b %d %Y %H %M %S")
    return datetime_a


def error_bar(ax, x, y, yerr, color, label=None):
    ax.plot(x, y, color, label=label)
    lower = [i - j for i, j in zip(y, yerr)]
    upper = [i + j for i, j in zip(y, yerr)]
    ax.fill_between(x, lower, upper, color=color, alpha=0.2, edgecolor='w')


def get_color(genotype):
    if "Cas3" in genotype:
        return "orangered"
    else:
        return "lightgrey"


def dot_plot(mice, func, ylabel):
    genotypes = {}
    for mouse in mice:
        if mouse.vector == "WT":
            continue
        hot_string = "%s+Hot" % mouse.vector
        cold_string = "%s+Cold" % mouse.vector
        if hot_string not in genotypes:
            genotypes[hot_string] = {}
        if mouse.mouse not in genotypes[hot_string]:
            genotypes[hot_string][mouse.mouse] = func(mouse.average_hot_pokes)
        if cold_string not in genotypes:
            genotypes[cold_string] = {}
        if mouse.mouse not in genotypes[cold_string]:
            genotypes[cold_string][mouse.mouse] = func(mouse.average_cold_pokes)
    labels = []
    last_genotype = None
    fig, ax = plt.subplots(nrows=2, gridspec_kw={"height_ratios": [1, 6]})
    fig.set_size_inches(6, 6)
    genotype_list = sorted(list(genotypes), reverse=True)
    num_statistically_significant = 0
    x_values = sorted(list(range(5, int(15 * len(genotype_list) / 2), 15)) + list(
        range(15, int(15 * len(genotype_list) / 2) + 1,
              15)))  # x values for where ticks go. (5, 15, 20, 30, 35, 45 ...)
    for i, genotype in enumerate(genotype_list):
        if last_genotype:
            if genotype.split("+")[0] != last_genotype.split("+")[0]:
                last_genotype = None
        pretty_genotype = genotype.replace("+", "\n").replace(" (", "\n").replace(")", "")
        if "hot" in genotype.lower():
            kwargs = {'facecolors': 'w'}
        else:
            kwargs = {}
        labels.append(pretty_genotype)
        total, count = 0, 0
        for mouse, value in genotypes[genotype].items():
            print(mouse, "%.2f" % value)
            if last_genotype is not None:
                if mouse in genotypes[last_genotype]:
                    last_value = genotypes[last_genotype][mouse]
                    ax[1].plot([x_values[i - 1], x_values[i]], [last_value, value], color=get_color(genotype),
                               zorder=2)
                    ax[1].scatter(x=x_values[i], y=value, color=get_color(genotype), s=50, zorder=3, **kwargs)
                    total += value
                    count += 1
            else:
                ax[1].scatter(x=x_values[i], y=value, color=get_color(genotype), s=50, zorder=3, **kwargs)
                total += value
                count += 1
        ax[1].hlines(total / count, x_values[i] - 2, x_values[i] + 2, zorder=1)
        last_genotype = genotype
        for j in range(i + 1, len(genotype_list)):
            first_arr = list(genotypes[genotype].values())
            second_arr = list(genotypes[genotype_list[j]].values())
            _, p = scipy.stats.ttest_ind(first_arr, second_arr, equal_var=False)
            if p < 0.2:
                ax[0].plot([x_values[i], x_values[i]],
                           [num_statistically_significant, num_statistically_significant + .5], "k-")
                ax[0].plot([x_values[j], x_values[j]],
                           [num_statistically_significant, num_statistically_significant + .5], "k-")
                ax[0].plot([x_values[i], x_values[j]],
                           [num_statistically_significant + .5, num_statistically_significant + .5], "k-")
                ax[0].text(x=(x_values[i] + x_values[j]) / 2, y=num_statistically_significant + .5, ha="center",
                           va="bottom", s="p=%.3f" % p, fontsize=12)
                num_statistically_significant += 2
        ax[0].set_ylim([0, num_statistically_significant + 1])

    ax[1].set_xticks(x_values)
    ax[1].set_xticklabels(labels)
    ax[1].set_ylabel(ylabel)
    ax[1].spines['top'].set_visible(False)
    ax[0].set_xlim(ax[1].get_xlim())
    ax[0].spines['bottom'].set_visible(False)
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    plt.subplots_adjust(hspace=0, wspace=0)
    plt.show()


if __name__ == "__main__":
    font_path = r"C:\Users\Fillan\Documents\Fonts\Myriad Pro Regular.ttf"
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = prop.get_name()
    plt.rcParams['font.size'] = 14
    plt.rcParams['axes.linewidth'] = 2

    HEIGHT = 1
    HEIGHT_BUFFER = 0.2
    keyfile = r"R:\Fillan\Parabrachial Ablations\Nose Poke\Keyfile.txt"
    folder = os.path.dirname(keyfile)
    mice = []
    with open(keyfile) as f:
        for line in f:
            mouse_number, vector = line.strip().split(":")
            mice.append(MouseData(mouse_number, vector))
    for mouse in mice:
        files = glob.glob(os.path.join(folder, "%s*.txt" % mouse.mouse))
        files = sorted(files, key=file_name_key)
        for i, file in enumerate(files):
            print(file)
            day_number = i // 2
            title = os.path.basename(file)[:-4]
            words = title.split(" ")
            dt = " ".join(words[1:-2])
            date = datetime.datetime.strptime(dt, "%b %d %Y %H %M").date()
            channel = words[-1]
            trial_start = None
            total_time_on = 0
            total_pokes = 0
            null = False
            with open(file) as f:
                for line in f:
                    if "NULL" in line:
                        print("NULL")
                        null = True
                        mouse.input(day_number, channel, None, None)
                        break
                    bout_start, bout_end = line.split(" - ")
                    bout_start = float(bout_start) / 60
                    bout_end = float(bout_end) / 60
                    if bout_end < 120:
                        total_time_on += (bout_end - bout_start)
                        total_pokes += 1
            if not null:
                mouse.input(day_number, channel, total_pokes, total_time_on)
        mouse.calculate_average_pokes()
    trial_labels = ["Total Pokes (Correct - Incorrect)", "Total Time On (Correct - Incorrect)"]
    trial_functions = [lambda x: x.pokes, lambda x: x.time_on]

    graph_titles = ["Cas3", "mCh", "WT"]
    for trial_label, trial_function in zip(trial_labels, trial_functions):
        fig, axs = plt.subplots(nrows=len(graph_titles), sharey=True)
        ymin, ymax = 0, 0
        for graph_title, ax in zip(graph_titles, axs):
            ax.set_xlabel("Day")
            ax.set_ylabel(trial_label)
            ax.set_title(graph_title)
            ax.plot([-1, TOTAL_DAYS+1], [0, 0], "k-", linewidth=3)
            ax.set_xlim([-0.5, TOTAL_DAYS - 0.5])
            for day_number in range(TOTAL_DAYS):
                active_values = []
                inactive_values = []
                for mouse in mice:
                    if mouse.vector == graph_title:
                        if not mouse.values[day_number].null:
                            active_values.append(trial_function(mouse.values[day_number])[0]-trial_function(mouse.values[day_number])[1])
                            inactive_values.append(0)
                ax.scatter([day_number] * len(active_values), active_values, facecolor='k', edgecolor='k')
                if len(active_values) > 0:
                    ax.bar([day_number], [sum(active_values) / len(active_values)], facecolor='none', edgecolor='k', width=1)
                ax.scatter([day_number] * len(inactive_values), inactive_values, facecolor='none', edgecolor='k')
                if len(inactive_values) > 0:
                    ax.bar([day_number], [sum(inactive_values) / len(inactive_values)], facecolor='none', edgecolor='k', width=1)
                this_y_min, this_y_max = ax.get_ylim()
                ymin = min(this_y_min, ymin)
                ymax = max(this_y_max, ymax)
        for j in range(len(graph_titles)):
            rect = Rectangle((-0.5, ymin), 1, ymax-ymin, edgecolor='none', facecolor=[1, 0, 0, 0.2])
            axs[j].add_patch(rect)
            rect = Rectangle((0.5, ymin), TOTAL_DAYS - 2, ymax-ymin, edgecolor='none', facecolor=[0, 0, 1, 0.2])
            axs[j].add_patch(rect)
            rect = Rectangle((TOTAL_DAYS - 1.5, ymin), 1, ymax-ymin, edgecolor='none', facecolor=[1, 0, 0, 0.2])
            axs[j].add_patch(rect)
            axs[j].set_ylim([ymin, ymax])
        plt.show()

    colors = ['orangered', 'lightgrey']
    fig, ax = plt.subplots()
    for graph_title, color in zip(graph_titles, colors):
        days, means, stds = [], [], []
        for day_number in range(TOTAL_DAYS - 2):
            active_values = []
            for mouse in mice:
                if mouse.vector == graph_title:
                    if not mouse.values[day_number].null:
                        active_values.append(mouse.values[day_number].pokes[0] - mouse.values[day_number].pokes[1])
            active_values = np.array(active_values)
            days.append(day_number)
            means.append(np.mean(active_values))
            stds.append(np.std(active_values))
        error_bar(ax, days, means, stds, color, graph_title)
    plt.xlabel("Day")
    plt.ylabel("Correct pokes - incorrect pokes")
    plt.legend()
    plt.show()

    dot_plot(mice, lambda x: x.pokes[0] - x.pokes[1], "Correct - incorrect pokes")
    dot_plot(mice, lambda x: x.pokes[0], "Correct pokes")
    dot_plot(mice, lambda x: x.time_on[0] - x.time_on[1], "Correct - incorrect time on")
    dot_plot(mice, lambda x: x.time_on[0], "Correct time on")

