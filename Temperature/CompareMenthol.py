import StressTest, DisplayLoggedTemperature
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import enum
import datetime
import time
from matplotlib import font_manager
from enum import Enum
import scipy.stats

sample_times = range(-60, 300, 15)

class Timer:
    def __init__(self):
        self.start_time = None
        self.elapsed = 0

    def start(self):
        assert self.start_time is None
        self.start_time = time.time()

    def stop(self):
        assert self.start_time is not None
        self.elapsed += time.time() - self.start_time
        self.start_time = None

    def __repr__(self):
        if self.start_time is None:
            return "%.2f" % self.elapsed
        else:
            return "%.2f" % self.elapsed + time.time() - self.start_time


class Mouse:
    def __init__(self, core_times, core_temperatures, baseline_temperature, mouse_number):
        self.core = []
        self.number = mouse_number
        for time in sample_times:
            self.core.append(Mouse.get_temperature_at_time(core_times, core_temperatures, time / 60))
        i = 0
        while core_times[i] < 0:
            i += 1
        self.average_core_temperature = sum(core_temperatures[i:]) / len(core_temperatures[i:]) + baseline_temperature

    @staticmethod
    def get_temperature_at_time(times, temperatures, time):
        assert len(times) == len(temperatures)
        if time < times[0]:
            return temperatures[0]
        if time > times[-1]:
            return temperatures[-1]
        for i in range(len(times) - 1):
            if times[i] < time <= times[i + 1]:
                return temperatures[i] + (temperatures[i + 1] - temperatures[i]) * (time - times[i])/ (times[i + 1] - times[i])


class MouseData:
    def __init__(self, line):
        self.number, self.vector = line.strip().split(":")
        self.drug_start = None
        self.control_start = None

    def check_date(self):
        if self.drug_start and self.control_start:
            assert self.drug_start.date() != self.control_start.date()

    def get_format(self):
        if self.vector == "Cas3":
            return "r"
        elif self.vector == "mCh":
            return "b"
        else:
            raise IOError("%s is not a valid genotype" % str)


def read_keyfile(key_file):
    all_mouse_numbers = {}
    with open(key_file) as f:
        currently_open_mouse_numbers = []
        reading_mouse_numbers = True
        for line in f:
            line = line.strip()
            if line[0] == ">":
                reading_mouse_numbers = False
                trial_type_str, trial_start = line[1:].split("-")
                if "menthol" in trial_type_str.lower() or "mk5046" in trial_type_str.lower():
                    for m in currently_open_mouse_numbers:
                        m.drug_start = datetime.datetime.strptime(trial_start, "%m/%d/%Y %I:%M %p")
                else:
                    for m in currently_open_mouse_numbers:
                        m.control_start = datetime.datetime.strptime(trial_start, "%m/%d/%Y %I:%M %p")
            elif line[0] == "#":
                continue
            else:
                if not reading_mouse_numbers:
                    for m in currently_open_mouse_numbers:
                        all_mouse_numbers[m.number] = m
                    currently_open_mouse_numbers = []
                    reading_mouse_numbers = True
                currently_open_mouse_numbers.append(MouseData(line))
        for m in currently_open_mouse_numbers:
            all_mouse_numbers[m.number] = m
    return all_mouse_numbers


def create_relative_temperatures(temperatures, start, relative_temp=False):
    import math
    delta_temperature_at_start = 0
    for dt, temp in temperatures:
        delta_temperature_at_start = temp
        if dt > start:
            break
    relative_times = []
    relative_temperatures = []
    for dt, temp in temperatures:
        if math.isnan(temp):
            continue
        if dt < start:
            relative_times.append(-1 * (start - dt).seconds / 3600)
        else:
            relative_times.append((dt - start).seconds / 3600)
        if relative_temp:
            relative_temperatures.append(temp - delta_temperature_at_start)
        else:
            relative_temperatures.append(temp)
    return relative_times, relative_temperatures, delta_temperature_at_start


def wobble_plot(dict):
    """
    :param dict: Ex: dict['Vglut2 Cas3'] = [(18, _), (19, _), (18.5, _)], dict['Vglut2 mCh'] = [(17, _)]
        where the values are the values at the start of the test
    :param ax: matplotlib axes
    :return:
    """
    fig, axs = plt.subplots(nrows=2, gridspec_kw={"height_ratios": [1, 6]})
    w = 0.8  # bar width
    sample_index = 10

    conditions = [(False, "mCh"), (False, "Cas3"), (False, ""), (True, "mCh"), (True, "Cas3")]
    num_statistically_significant = 0
    for x_position, condition in enumerate(conditions):
        received_drug, vector = condition
        if vector == "":
            continue
        ys = np.array([x.core[sample_index] for x in dict[received_drug][vector]])
        color = "orangered" if vector == "Cas3" else "grey"
        axs[1].bar(x_position,
               height=np.mean(ys), # error bars
               capsize=6,  # error bar cap width in points
               width=w,  # bar width,
               color=color,  # face color transparent
               edgecolor=(0, 0, 0, 0),
               # ecolor=colors,    # error bar colors; setting this raises an error for whatever reason.
               )
        axs[1].scatter([x_position] * len(ys), ys, edgecolor="k", facecolor=color)
        print(x_position, "%.2f±%.2f (n=%i)" % (np.mean(ys), np.std(ys), np.size(ys)))
        for new_x_position in range(x_position + 1, len(conditions)):
            new_received_drug, new_vector = conditions[new_x_position]
            if new_vector == "":
                continue
            new_ys = np.array([x.core[sample_index] for x in dict[new_received_drug][new_vector]])
            _, p = scipy.stats.ttest_ind(ys, new_ys, equal_var=False)
            if p < 1:
                axs[0].plot([x_position, x_position],
                            [num_statistically_significant, num_statistically_significant + .5],
                            "k-")
                axs[0].plot([new_x_position, new_x_position], [num_statistically_significant, num_statistically_significant + .5],
                            "k-")
                axs[0].plot([x_position, new_x_position],
                            [num_statistically_significant + .5, num_statistically_significant + .5], "k-")
                axs[0].text(x=(x_position + new_x_position) * 0.5, y=num_statistically_significant + .5, ha="center", va="bottom",
                            s="p=%.4f" % p, fontsize=12)
                num_statistically_significant += 2
    axs[0].set_ylim([0, num_statistically_significant + 1])
    axs[1].spines['top'].set_visible(False)
    axs[0].set_xlim(axs[1].get_xlim())
    axs[0].spines['bottom'].set_visible(False)
    axs[0].set_xticks([])
    axs[0].set_yticks([])
    plt.subplots_adjust(hspace=0, wspace=0)
    axs[1].axhline(y=0, xmin=0, xmax=6, color='k')
    axs[1].tick_params(axis='x', which='major', length=0)
    axs[1].set_xticks([0.5, 3.5])
    axs[1].set_xticklabels(["Vehicle", "Menthol"])
    axs[1].set_xticks([0, 1, 3, 4], minor=True)
    axs[1].set_xticklabels(["mCh", "Cas3", "mCh", "Cas3"], minor=True)
    axs[1].tick_params(axis='x', which='major', labelsize=16, pad=15)
    axs[1].tick_params(axis='x', which='minor', labelsize=12)
    axs[1].set_ylabel("ΔTemperature after %i minutes (°C)" % sample_times[sample_index])
    plt.show()


def dot_plot(dict, ax):
    """
    :param dict: Ex: dict['Vglut2 Cas3'] = [(18, 22), (19, 27), (18.5, 11)], dict['Vglut2 mCh'] = [(17, 17)]
        where the values are the values at the start of the test, then the value at the end of the test
    :param ax: matplotlib axes
    :return:
    """
    ticks = []
    for i, genotype in enumerate(list(dict.keys())):
        ticks.append(genotype.replace(" ", "\n"))
        for start, end in dict[genotype]:
            ax.scatter(i * 20, start, c='k')
            ax.scatter(i * 20 + 10, end, c='k')
            ax.plot([i * 20, i * 20 + 10], [start, end], MouseData.return_format(genotype))
    ax.set_xticks(range(5, len(ticks) * 20 + 5, 20))
    ax.set_xticklabels(ticks)
    ax.set_ylabel("Temperature")


def error_bar(ax, x, y, yerr, color, label=None):
    ax.plot(x, y, color, label=label)
    lower = [i - j for i, j in zip(y, yerr)]
    upper = [i + j for i, j in zip(y, yerr)]
    ax.fill_between(x, lower, upper, color=color, alpha=0.2, edgecolor='w')


def plot(mouse_data):
    temperatures = {}
    path = os.path.join(root_folder, "*", "*", "*.csv")
    for core_file in glob.glob(path):
        mouse_number = os.path.basename(core_file)[:-4]
        if mouse_number not in mouse_data:
            continue
        mouse_datum = mouse_data[mouse_number]
        mouse_datum.check_date()

        temps = DisplayLoggedTemperature.temperatures_from_file(core_file)
        core_temperatures = DisplayLoggedTemperature.average_temperatures(temps)
        date = core_temperatures[0][0].date()
        received_drug = False
        received_control = False
        if mouse_datum.drug_start:
            if mouse_datum.drug_start.date() == date:
                received_drug = True
        if mouse_datum.control_start:
            if mouse_datum.control_start.date() == date:
                received_control = True
        if received_drug and received_control:
            raise AssertionError("Mouse %s has two tests on %s" % (mouse_datum.number, mouse_datum.drug_start.date()))
        if not received_drug and not received_control:
            print("Mouse %s has no tests on %s" % (mouse_datum.number, date))
            continue
        core_times, core_temperatures, baseline_temp = create_relative_temperatures(core_temperatures, mouse_datum.drug_start if received_drug else mouse_datum.control_start, relative_temp=False)
        if received_drug not in temperatures:
            temperatures[received_drug] = {}
        if mouse_datum.vector not in temperatures[received_drug]:
            temperatures[received_drug][mouse_datum.vector] = []
        temperatures[received_drug][mouse_datum.vector].append(Mouse(core_times, core_temperatures, baseline_temp, mouse_number))
    fig, ax = plt.subplots()
    value_at_time = {"Cas3": [], "mCh": []}
    for vector in temperatures[False]:
        for mouse in temperatures[False][vector]:
            if vector == "Cas3":
                color = 'orangered'
            else:
                color = 'grey'
            menthol_mouse = None
            for i_mouse in temperatures[True][vector]:
                if mouse.number == i_mouse.number:
                    menthol_mouse = i_mouse
                    break
            if menthol_mouse is None:
                print("%s doesn't have a menthol trial" % mouse.number)
            core_values = []
            for i in range(len(sample_times)):
                core_values.append(menthol_mouse.core[i] - mouse.core[i])
                if i == 9:
                    value_at_time[vector].append(menthol_mouse.core[i] - mouse.core[i])
            ax.plot(sample_times, core_values, color=color)
    plt.xlim([-60, 180])
    ax.axvline(x=0, c='k')
    ax.set_ylabel("ΔCore Temp (C)")
    ax.set_xticks([-60, 0, 60, 120, 180])
    ax.set_xlabel("Minutes from trial start")
    plt.title("Menthol - Vehicle")
    plt.legend(prop={"size": 10})
    plt.show()

    for i, time in enumerate(sample_times):
        print(i, time)

    plt.scatter([10] * len(value_at_time["Cas3"]), value_at_time["Cas3"], color="orangered")
    plt.scatter([20] * len(value_at_time["mCh"]), value_at_time["mCh"], color="gray")
    _, p = scipy.stats.ttest_ind(value_at_time["Cas3"], value_at_time["mCh"])
    plt.title("p=%.4f" % p)
    plt.show()

    for received_drug, tests in temperatures.items():
        fig, ax = plt.subplots()
        for vector, mice in tests.items():
            core_means, core_stds = [], []
            for i in range(len(sample_times)):
                if vector == "Cas3":
                    color = 'orangered'
                else:
                    color = 'grey'
                core_values = []
                for mouse in mice:
                    core_values.append(mouse.core[i])
                    #ax.scatter(sample_times[i], mouse.core[i], color=color)
                core_values = np.array(core_values)
                core_means.append(np.mean(core_values))
                core_stds.append(np.std(core_values))
            error_bar(ax, x=sample_times, y=core_means, yerr=core_stds, color=color, label="%s (n=%s)" % (vector, len(mice)))
        plt.xlim([-60, 180])
        ax.axvline(x=0, c='k')
        ax.set_ylabel("ΔCore Temp (C)")
        ax.set_xticks([-60, 0, 60, 120, 180])
        ax.set_xlabel("Minutes from trial start")
        plt.title("Menthol" if received_drug else "Vehicle")
        plt.legend(prop={"size":10})
        plt.show()
    wobble_plot(temperatures)


if __name__ == "__main__":
    font_path = r"C:\Users\Fillan\Documents\Fonts\Myriad Pro Regular.ttf"
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = prop.get_name()
    plt.rcParams['font.size'] = 14
    plt.rcParams['axes.linewidth'] = 2

    root_folder = r"R:\Fillan\Parabrachial Ablations\Menthol"
    key_file = os.path.join(root_folder, "Keyfile.txt")
    mouse_numbers = read_keyfile(key_file)
    plot(mouse_numbers)