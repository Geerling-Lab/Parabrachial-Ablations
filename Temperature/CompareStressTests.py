import StressTest, DisplayLoggedTemperature
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import enum
import datetime
import time
from matplotlib import font_manager
import scipy.stats

sample_times = range(-60, 315, 15)

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


class TrialType(enum.Enum):
    Control = 1
    Hot = 2
    Cold = 3

    @staticmethod
    def get_type(str):
        str = str.lower()
        if "control" in str:
            return TrialType.Control
        elif "hot" in str:
            return TrialType.Hot
        elif "cold" in str:
            return TrialType.Cold
        else:
            raise IOError("Invalid trial type: %s" % str)

    @staticmethod
    def get_str(trial_type):
        if trial_type == TrialType.Control:
            return "Control"
        elif trial_type == TrialType.Hot:
            return "Hot"
        elif trial_type == TrialType.Cold:
            return "Cold"


class MouseData:
    def __init__(self, line):
        self.number, self.genotype, self.vector = line.strip().split(":")
        self.start = {}

    def get_format(self):
        return MouseData.return_format("%s %s" % (self.genotype, self.vector))

    def get_color(self):
        return MouseData.return_color(self.genotype)

    def get_linestyle(self):
        if "mCh" in self.vector:
            return "--"
        elif "Cas3" in self.vector:
            return "-"
        else:
            raise IOError("%s is not a valid vector" % self.vector)

    @staticmethod
    def return_color(str):
        if str == "Vglut2":
            return "r"
        elif str == "FoxP2":
            return "orangered"
        elif str == "Lmx1b":
            return "violet"
        elif str == "Pdyn":
            return "b"
        elif str == "Sst":
            return "c"
        elif str == "Cck":
            return "m"
        elif str == "NPS":
            return "darkorange"
        else:
            raise IOError("%s is not a valid genotype" % str)


class Mouse:
    def __init__(self, env_times, env_temperatures, core_times, core_temperatures, baseline_temperature):
        self.environment, self.core = [], []
        for time in sample_times:
            self.environment.append(Mouse.get_temperature_at_time(env_times, env_temperatures, time / 60))
            self.core.append(Mouse.get_temperature_at_time(core_times, core_temperatures, time / 60))
        i = 0
        while core_times[i] < 0:
            i += 1
        self.average_core_temperature = sum(core_temperatures[i:]) / len(core_temperatures[i:]) + baseline_temperature
        self.baseline_temperature = baseline_temperature

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


def read_keyfile(key_file):
    all_mouse_numbers = {}
    with open(key_file) as f:
        currently_open_mouse_numbers = []
        reading_mouse_numbers = True
        for line in f:
            line = line.strip()
            if line[0] == ">":
                reading_mouse_numbers = False
                trial_type_str, trial_start = line[1:].split(" ")
                trial_type = TrialType.get_type(trial_type_str)
                for m in currently_open_mouse_numbers:
                    m.start[trial_type] = datetime.datetime.strptime(trial_start, "%H:%M")
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


def wobble_plot(dict, ax, func):
    """
    :param dict: Ex: dict['Vglut2] = {'Cas3': Mouse1, Mouse2, Mouse3}, {'mCh': Mouse4}
        where the values are the values at the start of the test
    :param ax: matplotlib axes
    :return:
    """
    w = 0.8  # bar width
    labels = [("Vglut2","mCh"), ("Vglut2","Cas3"), "",
              ("FoxP2","mCh"), ("FoxP2","Cas3"), "",
              ("Lmx1b","mCh"), ("Lmx1b","Cas3"), "",
              ("Pdyn","mCh"), ("Pdyn","Cas3"), "",
              ("Sst","mCh"), ("Sst","Cas3"), "",
              ("Cck","mCh"), ("Cck","Cas3"), "",
              ("NPS","mCh"), ("NPS", "Cas3"), ""]
    vector_labels = []
    vector_labels_x = []
    for i, label in enumerate(labels):
        if len(label) > 1:
            vector_labels.append(label[1])
            vector_labels_x.append(i)
    x = list(range(len(labels)))
    y = []
    for l in labels:
        if len(l) < 1:
            y.append([])
        else:
            if l[0] in dict:
                if l[1] in dict[l[0]]:
                    y.append([func(i) for i in dict[l[0]][l[1]]])
                else:
                    y.append([])
            else:
                y.append([])

    colors = [(1, 0, 0, 0.2) if "Cas3" in x else (.6, .6, .6, 0.4) for x in labels]
    ax.bar(x,
           height=[np.mean(yi) for yi in y],
           width=w,  # bar width,
           color=colors,  # face color transparent
           edgecolor=(0, 0, 0, 0),
           # ecolor=colors,    # error bar colors; setting this raises an error for whatever reason.
           )
    ax.axhline(y=0, xmin=0, xmax=1, color='k')

    for i, label in enumerate(labels):
        # distribute scatter randomly across whole width of bar
        color = (1, 0, 0) if "Cas3" in label else (0.6, 0.6, 0.6)
        ax.scatter(x[i] + np.random.random(len(y[i])) * w / 8 - w / 16, y[i], edgecolor="k", facecolor=color)
        print("%s: Mean=%.2f±%.2f (n=%i)" % (label, np.mean(y[i]), np.std(y[i]), len(y[i])))
        for j in range(i+1, len(labels)):
            if len(labels[i]) > 0 and len(labels[j]) > 0:
                if labels[i][0] == labels[j][0]:
                    _, p = scipy.stats.ttest_ind(y[i], y[j])
                    print("%s: p=%.4f" % (labels[i][0], p))

    ax.tick_params(axis='x', which='major', length=0)
    ax.set_xticks([l + 0.5 for l in x[::3]])
    ax.set_xticklabels([l[0] for l in labels[::3]])
    ax.set_xticks(vector_labels_x, minor=True)
    ax.set_xticklabels(vector_labels, minor=True)
    ax.tick_params(axis='x', which='major', labelsize=14, pad=15)
    ax.tick_params(axis='x', which='minor', labelsize=8)


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


def scatter_plot(dict, time, str_time):
    fig, ax = plt.subplots(nrows=2)
    control_cores, control_environments = [], []
    for genotype, vector_dict in dict.items():
        for vector, mice in vector_dict.items():
            color = 'w'
            if "Cas3" in vector:
                color = MouseData.return_color(genotype)
            cores, environments = [], []
            for mouse in mice:
                cores.append(mouse.core[time] + mouse.baseline_temperature)
                environments.append(mouse.environment[time])
                if "mCh" in vector:
                    control_cores.append(cores[-1])
                    control_environments.append(environments[-1])
            ax[0].scatter(environments, cores, color=color, edgecolors=MouseData.return_color(genotype))
    control_cores = np.array(control_cores)
    control_environments = np.array(control_environments)
    m, b, r, *_ = scipy.stats.linregress(control_environments, control_cores)
    _, p = scipy.stats.pearsonr(control_environments, control_cores)
    ax[0].plot(control_environments, np.poly1d([m, b])(control_environments), 'k', label="p=%.5f" % p)
    ax[0].legend()
    ax[0].set_ylabel("Core Temperature after %s" % str_time)
    ax[0].set_xlabel("Environment Temperature after %s" % str_time)
    ax[1].set_ylabel("ΔCore Temperature from trendline")
    func = lambda x: x.core[time] - m * x.environment[time] - b + x.baseline_temperature
    wobble_plot(dict, ax[1], func)
    plt.show()


def plot_condition(mouse_data, trial_type):
    temperatures = {}
    f = open(r"R:\Fillan\Parabrachial Ablations\Temperature\Temperature at Times.csv", 'w')
    for core_file in glob.glob(os.path.join(root_folder, "*", "%s *" % TrialType.get_str(trial_type), "Core *.csv")):
        folder = os.path.dirname(core_file)
        mouse_number = os.path.basename(core_file).split(" ")[1][:-4]
        if "repeat" in core_file.lower():
            continue
        mouse_datum = mouse_data[mouse_number]
        env_file = os.path.join(folder, "Elitech %s.txt" % mouse_number)
        if not os.path.exists(env_file):
            print("Couldn't find environmental file for %s" % core_file)
            continue
        env_temperatures = StressTest.read_elitech_file(env_file)
        if trial_type == TrialType.Hot:
            for dt, temp in env_temperatures:
                if temp > 34:
                    break
            else:
                print("%s didn't get hot enough" % mouse_number)
                #continue
        temps = DisplayLoggedTemperature.temperatures_from_file(core_file)
        core_temperatures = DisplayLoggedTemperature.average_temperatures(temps)
        date = core_temperatures[0][0].date()
        start_time = datetime.datetime.combine(date, mouse_datum.start[trial_type].time())
        core_times, core_temperatures, baseline_temp = create_relative_temperatures(core_temperatures, start_time, relative_temp=True)
        env_times, env_temperatures, _ = create_relative_temperatures(env_temperatures, start_time, relative_temp=False)
        if mouse_datum.genotype not in temperatures:
            temperatures[mouse_datum.genotype] = {}
        if mouse_datum.vector not in temperatures[mouse_datum.genotype]:
            temperatures[mouse_datum.genotype][mouse_datum.vector] = []
        temperatures[mouse_datum.genotype][mouse_datum.vector].append(Mouse(env_times, env_temperatures, core_times, core_temperatures, baseline_temp))

    slopes = []
    intercepts = []
    for i in range(len(sample_times)):
        core_values = []
        env_values = []
        for genotype, vectors in temperatures.items():
            for vector, mice in vectors.items():
                if "mCh" in vector:
                    for mouse in mice:
                        core_values.append(mouse.core[i] + mouse.baseline_temperature)
                        env_values.append(mouse.environment[i])
        m, b, r, *_ = scipy.stats.linregress(env_values, core_values)
        slopes.append(m)
        intercepts.append(b)

    for genotype, vectors in temperatures.items():
        fig, ax = plt.subplots(nrows=2, sharex=True, gridspec_kw={'height_ratios': [1, 2]})
        plt.subplots_adjust(wspace=0, hspace=0)
        for vector, mice in vectors.items():
            core_means, core_stds, env_means, env_stds = [], [], [], []
            for i in range(len(sample_times)):
                core_values, env_values = [], []
                for mouse in mice:
                    core_values.append(mouse.core[i] + mouse.baseline_temperature - mouse.environment[i] * slopes[i] - intercepts[i])
                    env_values.append(mouse.environment[i])
                core_values = np.array(core_values)
                env_values = np.array(env_values)
                core_means.append(np.mean(core_values))
                core_stds.append(np.std(core_values))
                env_means.append(np.mean(env_values))
                env_stds.append(np.std(env_values))
            if "mch" in vector.lower():
                color = 'slategray'
                fmt = "b-"
            else:
                color = MouseData.return_color(genotype)
                fmt = "r-"
            error_bar(ax[0], x=sample_times, y=env_means, yerr=env_stds, color=color, label="%s (n=%s)" % (vector, len(mice)))
            error_bar(ax[1], x=sample_times, y=core_means, yerr=core_stds, color=color)
            #ax[0].errorbar(x=sample_times, y=env_means, yerr=env_stds, fmt=fmt)
            #ax[1].errorbar(x=sample_times, y=core_means, yerr=core_stds, fmt=fmt)
        ax[1].set_xlim([-60, 300])
        ax[1].set_ylim([-2, 2])
        ax[1].axvline(x=0, c='k')
        ax[1].set_ylabel("ΔCore Temp (C)")
        ax[1].set_xticks([-60, 0, 60, 120, 180, 240])
        ax[0].set_xlim([-60, 300])
        ax[0].axvline(x=0, c='k')
        ax[0].set_yticks([5, 22])
        ax[0].set_ylabel("Cage Temp (C)")
        ax[1].set_xlabel("Minutes from trial start")
        fig.suptitle(genotype)
        ax[0].legend(prop={"size":10})
        plt.show()

    time, str_time = 20, "4 hours"
    if trial_type == TrialType.Hot:
        time, str_time = 10, "1.5 hour"

    scatter_plot(temperatures, time=time, str_time=str_time)

    fig, ax = plt.subplots()
    print("Average core temperature during test")
    wobble_plot(temperatures, ax, lambda x: x.average_core_temperature)
    ax.set_ylabel("Average core temperature during test")
    plt.show()
    fig, ax = plt.subplots()
    print("ΔCore temperature (C) over 4 hours")
    wobble_plot(temperatures, ax, lambda x: x.core[20])
    ax.set_ylabel("ΔCore temperature (C) over 4 hours")
    plt.show()
    fig, ax = plt.subplots()
    print("ΔCore temperature (C) over 60 minutes")
    wobble_plot(temperatures, ax, lambda x: x.core[8])
    ax.set_ylabel("ΔCore temperature (C) over 60 minutes")
    plt.show()


if __name__ == "__main__":
    font_path = r"C:\Users\Fillan\Documents\Fonts\Myriad Pro Regular.ttf"
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = prop.get_name()
    plt.rcParams['font.size'] = 14
    plt.rcParams['axes.linewidth'] = 2

    root_folder = r"R:\Fillan\Parabrachial Ablations\Temperature"
    key_file = os.path.join(root_folder, "Keyfile.txt")
    mouse_numbers = read_keyfile(key_file)
    plot_condition(mouse_numbers, TrialType.Hot)
    plot_condition(mouse_numbers, TrialType.Cold)
    plot_condition(mouse_numbers, TrialType.Control)