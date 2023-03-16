import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from enum import Enum
import datetime
import math
from matplotlib import font_manager
from bee_swarm import fillan_bee_swarm
import scipy.stats


class Vector(Enum):
    Cas3 = 1
    mCh = 2


class Genotype(Enum):
    Vglut2 = 1
    Foxp2 = 2
    Lmx1b = 3
    Pdyn = 4
    Cck = 5
    Sst = 6
    NPS = 7
    Unknown = 8


class Mouse:
    def __init__(self, number, genotype, vector, birth_date, dates, weights, male_sex, stereology=None):
        self.genotype = genotype
        self.vector = vector
        self.number = number
        self.birth_date = birth_date
        self.male_sex = male_sex
        self.dates = []
        self.weights = []
        self.stereology = None
        if stereology:
            if type(stereology) == str:
                self.stereology = None
            else:
                if not math.isnan(stereology):
                    self.stereology = stereology
        surgery_weight = weights[0]
        surgery_day = dates[0]
        for weight in weights:
            if not math.isnan(weight):
                self.weights.append(weight - surgery_weight)
        for date in dates:
            delta_days = (date - surgery_day).days
            if delta_days < 0:
                print(self.genotype, self.number)
            assert type(delta_days) == int
            self.dates.append(delta_days)
        if len(self.weights) != len(self.dates):
            raise IOError()

    def draw(self, ax):
        if self.vector == Vector.Cas3:
            color = 'r'
        else:
            color = 'grey'
        ax.plot(self.dates, self.weights, color)

    def get_weight_change(self):
        if len(self.weights) < 2:
            return 0
        return self.weights[-1] - self.weights[1]

    def get_weight_at_date(self, date):
        smoothing = 15
        for i in range(len(self.dates)):
            if self.dates[i] > date:
                if self.dates[i - 1] == self.dates[i]:
                    print("Two weights on same date", mouse_number, genotype, len(self.dates))
                slope = (self.weights[i] - self.weights[i - 1]) / (self.dates[i] - self.dates[i - 1])
                weight = slope * (date - self.dates[i - 1]) + self.weights[i - 1]
                certainty = 1 / (min(date - self.dates[i - 1], self.dates[i] - date) + smoothing)
                break
        else:
            weight = self.weights[-1]
            certainty = 1 / ((date - self.dates[-1]) + smoothing)
        certainty = 1
        return weight, certainty


def plot_individual_mice(mice, func, title):
    fig, ax = plt.subplots()
    for mouse in mice:
        if func(mouse):
            mouse.draw(ax)
    plt.xlabel("Days from birth")
    plt.ylabel("Weight")
    plt.xlim([0, 100])
    plt.title(title)
    plt.show()


def error_bar(ax, x, y, yerr, color, label=None):
    ax.plot(x, y, color, label=label)
    lower = [i - j for i, j in zip(y, yerr)]
    upper = [i + j for i, j in zip(y, yerr)]
    ax.fill_between(x, lower, upper, color=color, alpha=0.2, edgecolor='w')


def plot_averaged_lines(mice, func, title, color='r'):
    fig, ax = plt.subplots()
    Cas3_n = 0
    mCh_n = 0
    xs, mCh_ys, Cas3_ys, mCh_yerr, Cas3_yerr = [], [], [], [], []
    for x in range(1, 200):
        xs.append(x)
        mCh_weights, mCh_certainties = [], []
        Cas3_weights, Cas3_certainties = [], []
        for mouse in mice:
            if func(mouse):
                if len(mouse.weights) >= 2:
                    weight, certainty = mouse.get_weight_at_date(x)
                    if mouse.vector == Vector.mCh:
                        mCh_weights.append(weight)
                        mCh_certainties.append(certainty)
                    else:
                        Cas3_weights.append(weight)
                        Cas3_certainties.append(certainty)
        Cas3_n = len(Cas3_weights)
        mCh_n = len(mCh_weights)
        mCh_ys.append(sum([weight * certainty for weight, certainty in zip(mCh_weights, mCh_certainties)]) / sum(mCh_certainties))
        mCh_yerr.append(np.sqrt(np.cov(mCh_weights, aweights=mCh_certainties)))
        Cas3_ys.append(sum([weight * certainty for weight, certainty in zip(Cas3_weights, Cas3_certainties)]) / sum(Cas3_certainties))
        Cas3_yerr.append(np.sqrt(np.cov(Cas3_weights, aweights=Cas3_certainties)))
    error_bar(ax, x=xs, y=mCh_ys, yerr=mCh_yerr, color='gray', label="n=%i" % mCh_n)
    error_bar(ax, x=xs, y=Cas3_ys, yerr=Cas3_yerr, color=color, label="n=%i" % Cas3_n)
    for mouse in mice:
        if func(mouse):
            if mouse.vector == Vector.mCh:
                #ax.plot(mouse.dates[1], mCh_ys[mouse.dates[1]], color="gray", marker="o", ms=6, zorder=1)
                if mouse.dates[-1] < len(mCh_ys):
                    ax.scatter(mouse.dates[-1], mCh_ys[mouse.dates[-1]], color="gray", marker="+", s=50, linewidth=2, zorder=1)
            else:
                #ax.plot(mouse.dates[1], Cas3_ys[mouse.dates[1]], color=color, marker="o", ms=6, zorder=1)
                if mouse.dates[-1] < len(Cas3_ys):
                    ax.scatter(mouse.dates[-1], Cas3_ys[mouse.dates[-1]], color=color, marker="+", s=50, linewidth=2, zorder=1)
    plt.xlabel("Days after surgery")
    plt.ylabel("Weight change")
    plt.xlim([0, 100])
    plt.ylim([-2, 8])
    plt.legend()
    plt.title(title)
    plt.hlines(0, 0, 100, color='k')
    plt.show()


def plot_bar_graph(mice, func, title, color='r'):
    fig, ax = plt.subplots(nrows=2, gridspec_kw={"height_ratios": [1, 6]})
    fig.set_size_inches(6, 8)
    Cas3, mCh = [], []
    for mouse in mice:
        if func(mouse):
            if mouse.dates[-1] - mouse.dates[0] > 30:
                value = (mouse.weights[-1] - mouse.weights[0]) / (mouse.dates[-1] - mouse.dates[0])
                print(value)
                if mouse.vector == Vector.Cas3:
                    Cas3.append(value)
                else:
                    mCh.append(value)
    ax[1].scatter(fillan_bee_swarm(mCh, 1, 0.5, 0.1, None), mCh, color='gray', zorder=2)
    ax[1].scatter(fillan_bee_swarm(Cas3, 3, 0.5, 0.1, None), Cas3, color=color, zorder=2)
    ax[1].plot([.5, 1.5], [sum(mCh) / len(mCh), sum(mCh) / len(mCh)], "k-", zorder=1)
    ax[1].plot([2.5, 3.5], [sum(Cas3) / len(Cas3), sum(Cas3) / len(Cas3)], "k-", zorder=1)
    ax[1].set_ylabel("Post-surgery weight change (g/day)")
    ax[1].set_xticks([1, 3], ["%s\nmCh" % title, "%s\nCas3" % title])
    _, p = scipy.stats.ttest_ind(Cas3, mCh)
    mCh = np.array(mCh)
    Cas3 = np.array(Cas3)
    print("mCh: %.3f ± %.3f, Cas3: %.3f ± %.3f" % (np.mean(mCh), np.std(mCh), np.mean(Cas3), np.std(Cas3)))
    print("mCh n=%i, Casp3 n=%i" % (len(mCh), len(Cas3)))
    ax[0].plot([1, 1], [0, 0.5], "k-")
    ax[0].plot([3, 3], [0, 0.5], "k-")
    ax[0].plot([1, 3], [0.5, 0.5], "k-")
    ax[0].text(x=2, y=0.5, ha="center", va="bottom", s="p=%.4f" % p, fontsize=12)
    xmin, xmax = ax[1].get_xlim()
    ax[1].plot([0, 10], [0, 0], 'k-', zorder=0, linewidth=2)
    ax[1].set_xlim([xmin, xmax])
    ax[1].set_ylim([-.2, .2])
    ax[0].set_xticks([])
    ax[0].set_yticks([])
    ax[0].set_ylim([0, 2])
    ax[1].spines['top'].set_visible(False)
    ax[0].set_xlim(ax[1].get_xlim())
    ax[0].spines['bottom'].set_visible(False)
    plt.subplots_adjust(hspace=0, wspace=0)
    plt.show()


if __name__ == "__main__":
    font_path = r"C:\Users\Fillan\Documents\Fonts\Myriad Pro Regular.ttf"
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = prop.get_name()
    plt.rcParams['font.size'] = 18
    plt.rcParams['axes.linewidth'] = 2

    files = {r"R:\Fillan\Parabrachial Ablations\Lesion Mice Overview.xlsx": Genotype.Vglut2,
             r"R:\Fillan\Parabrachial Ablations\Experiment 3 Vglut2 with EEG\Vglut2 with EEG mice.xlsx": Genotype.Vglut2,
             r"R:\Fillan\Parabrachial Ablations\Experiment 3 Vglut2 without EEG\Vglut2 without EEG mice.xlsx": Genotype.Vglut2,
             r"R:\Fillan\Parabrachial Ablations\Experiment 3 Foxp2\FoxP2 mice.xlsx": Genotype.Foxp2,
             r"R:\Fillan\Parabrachial Ablations\Experiment 3 Lmx1b\Lmx1b mice.xlsx": Genotype.Lmx1b,
             r"R:\Fillan\Parabrachial Ablations\Subtypes of Foxp2 for Temperature.xlsx": Genotype.Unknown}
    mice = []
    for file_name, genotype in files.items():
        print(file_name)
        overview_df = pd.read_excel(file_name, sheet_name="Overview")
        weight_df = pd.read_excel(file_name, sheet_name="Weight", date_parser=['date'])
        mice_numbers = list(weight_df.columns)
        for index in range(0, len(mice_numbers), 2):
            mouse_number = mice_numbers[index]
            dates = [x.date() for x in list(weight_df[mice_numbers[index]]) if x is not pd.NaT]
            weights = [x for x in list(weight_df[mice_numbers[index + 1]]) if not math.isnan(x)]
            assert len(dates) == len(weights)

            try:
                row = overview_df.loc[overview_df["Mouse"] == int(mouse_number)]
            except ValueError as ve:
                row = overview_df.loc[overview_df["Mouse"] == mouse_number]
            if len(row) == 0:
                raise IOError("%s not found in overview" % mouse_number)
            if "Cas3" in row["Vector"].values[0]:
                vector = Vector.Cas3
            elif "mCh" in row["Vector"].values[0]:
                vector = Vector.mCh
            else:
                raise IOError("%s is not a valid vector" % row["Vector"])
            if genotype == Genotype.Unknown:
                if "Cck" in row["Genotype"].values[0]:
                    new_genotype = Genotype.Cck
                elif "Pdyn" in row["Genotype"].values[0]:
                    new_genotype = Genotype.Pdyn
                elif "NPS" in row["Genotype"].values[0]:
                    new_genotype = Genotype.NPS
                elif "Sst" in row["Genotype"].values[0]:
                    new_genotype = Genotype.Sst
                else:
                    raise IOError("%s is not a valid genotype" % row["Genotype"])
            else:
                new_genotype = genotype
            male_sex = True
            if "F" == row["Sex"].values[0]:
                male_sex = False
            if "Stereology" in row:
                stereology = row["Stereology"].values[0]
            else:
                stereology = None
            birth_date = pd.Timestamp(row["Date of Birth"].values[0]).to_pydatetime().date()
            mice.append(Mouse(mouse_number, new_genotype, vector, birth_date, dates, weights, male_sex, stereology))

    plot_bar_graph(mice, lambda x: x.genotype == Genotype.Vglut2, title="Vglut2")
    plot_individual_mice(mice, lambda x: x.genotype == Genotype.Vglut2, title="Vglut2")
    plot_averaged_lines(mice, lambda x: x.genotype == Genotype.Vglut2, title="Vglut2")
    plot_bar_graph(mice, lambda x: x.genotype == Genotype.Foxp2, title="FoxP2", color="orangered")
    plot_individual_mice(mice, lambda x: x.genotype == Genotype.Foxp2, title="Foxp2")
    plot_averaged_lines(mice, lambda x: x.genotype == Genotype.Foxp2, title="Foxp2", color="orangered")
    plot_bar_graph(mice, lambda x: x.genotype == Genotype.Lmx1b, title="Lmx1b", color="violet")
    plot_individual_mice(mice, lambda x: x.genotype == Genotype.Lmx1b, title="Lmx1b")
    plot_averaged_lines(mice, lambda x: x.genotype == Genotype.Lmx1b, title="Lmx1b", color="violet")
    plot_bar_graph(mice, lambda x: x.genotype == Genotype.Cck, title="Cck")
    plot_individual_mice(mice, lambda x: x.genotype == Genotype.Cck, title="Cck")
    plot_averaged_lines(mice, lambda x: x.genotype == Genotype.Cck, title="Cck")
    plot_bar_graph(mice, lambda x: x.genotype == Genotype.Pdyn, title="Pdyn")
    plot_individual_mice(mice, lambda x: x.genotype == Genotype.Pdyn, title="Pdyn")
    plot_averaged_lines(mice, lambda x: x.genotype == Genotype.Pdyn, title="Pdyn")
    plot_bar_graph(mice, lambda x: x.genotype == Genotype.NPS, title="NPS")
    plot_individual_mice(mice, lambda x: x.genotype == Genotype.NPS, title="NPS")
    plot_averaged_lines(mice, lambda x: x.genotype == Genotype.NPS, title="NPS")
    plot_bar_graph(mice, lambda x: x.genotype == Genotype.Sst, title="Sst")
    plot_individual_mice(mice, lambda x: x.genotype == Genotype.Sst, title="Sst")
    plot_averaged_lines(mice, lambda x: x.genotype == Genotype.Sst, title="Sst")

    stereological_counts = []
    weight_changes = []
    colors = []
    for mouse in mice:
        if mouse.genotype == Genotype.Vglut2:
            if mouse.stereology:
                stereological_counts.append(mouse.stereology)
                weight_changes.append(mouse.get_weight_change())
                if mouse.vector == Vector.Cas3:
                    colors.append('r')
                else:
                    colors.append('b')
    plt.scatter(stereological_counts, weight_changes, color=colors)
    plt.plot(np.unique(stereological_counts), np.poly1d(np.polyfit(stereological_counts, weight_changes, 1))(np.unique(stereological_counts)), 'k-')
    plt.xlabel("Stereological Count")
    plt.ylabel("Weight Change (g)")
    plt.tight_layout()
    plt.show()