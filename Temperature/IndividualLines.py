import StressTest, DisplayLoggedTemperature
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import enum
import datetime
import time
from matplotlib import font_manager


class Line_Event_Handler:
    def __init__(self, fig, ax):
        self.lines = []
        self.names = []
        self.ax = ax
        self.fig = fig
        self.annot = self.ax.annotate("", xy=(0, 0), xytext=(-20, 20), textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="w"),
                            arrowprops=dict(arrowstyle="->"))
        self.annot.set_visible(False)

    def add_line(self, line):
        self.lines.append(line)

    def hover(self, event):
        vis = self.annot.get_visible()
        if event.inaxes == self.ax:
            for line in self.lines:
                cont, ind = line.contains(event)
                if cont:
                    self.annot.xy = event.xdata, event.ydata
                    self.annot.set_text(line.get_label())
                    self.annot.get_bbox_patch().set_alpha(0.4)
                    self.annot.set_visible(True)
                    self.fig.canvas.draw_idle()
                else:
                    if vis:
                        self.annot.set_visible(False)
                        self.fig.canvas.draw_idle()


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
            return "g"
        elif str == "FoxP2":
            return "r"
        elif str == "Lmx1b":
            return "y"
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
    :param dict: Ex: dict['Vglut2 Cas3'] = [(18, _), (19, _), (18.5, _)], dict['Vglut2 mCh'] = [(17, _)]
        where the values are the values at the start of the test
    :param ax: matplotlib axes
    :return:
    """
    w = 0.8  # bar width
    labels = ["Vglut2 mCh", "Vglut2 Cas3", "",
              "FoxP2 mCh", "FoxP2 Cas3", "",
              "Lmx1b mCh", "Lmx1b Cas3", "",
              "Pdyn mCh", "Pdyn Cas3", "",
              "Sst mCh", "Sst Cas3", "",
              "Cck mCh", "Cck Cas3", "",
              "NPS mCh", "NPS Cas3", ""]
    vector_labels = []
    vector_labels_x = []
    for i, label in enumerate(labels):
        if " " in label:
            vector_labels.append(label.split(" ")[1])
            vector_labels_x.append(i)
    x = list(range(len(labels)))
    y = []
    for l in labels:
        if l not in dict:
            y.append([])
        else:
            y.append([func(i) for i in dict[l]])

    colors = [(1, 0, 0, 0.2) if "Cas3" in x else (0, 0, 1, 0.2) for x in labels]
    ax.bar(x,
           height=[np.mean(yi) for yi in y],
           yerr=[np.std(yi) for yi in y],  # error bars
           capsize=6,  # error bar cap width in points
           width=w,  # bar width,
           color=colors,  # face color transparent
           edgecolor=(0, 0, 0, 0),
           # ecolor=colors,    # error bar colors; setting this raises an error for whatever reason.
           )
    ax.axhline(y=0, xmin=0, xmax=1, color='k')

    for i, label in enumerate(labels):
        # distribute scatter randomly across whole width of bar
        color = (1, 0, 0) if "Cas3" in label else (0, 0, 1)
        ax.scatter(x[i] + np.random.random(len(y[i])) * w / 8 - w / 16, y[i], edgecolor="k", facecolor=color)

    ax.tick_params(axis='x', which='major', length=0)
    ax.set_xticks([l + 0.5 for l in x[::3]])
    ax.set_xticklabels([l.split(" ")[0] for l in labels[::3]])
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
            ax.plot([i * 20, i * 20 + 10], [start, end], c=MouseData.return_color(genotype))
    ax.set_xticks(range(5, len(ticks) * 20 + 5, 20))
    ax.set_xticklabels(ticks)
    ax.set_ylabel("Temperature")


def plot_condition(mouse_data, trial_type):
    fig, ax = plt.subplots(nrows=2, sharex=True, gridspec_kw={'height_ratios': [1, 2]})
    plt.subplots_adjust(wspace=0, hspace=0)
    line_event_handlers = np.empty(ax.shape, dtype=object)
    for i in range(2):
        line_event_handlers[i] = Line_Event_Handler(fig, ax[i])
    baseline_temps = {}
    f = open(r"R:\Fillan\Parabrachial Ablations\Temperature\Temperature at Times.csv", 'w')
    for core_file in glob.glob(os.path.join(root_folder, "*", "%s *" % TrialType.get_str(trial_type), "Core *.csv")):
        folder = os.path.dirname(core_file)
        mouse_number = os.path.basename(core_file).split(" ")[1][:-4]
        if "repeat" in core_file.lower():
            continue
        mouse_datum = mouse_data[mouse_number]
        if mouse_datum.genotype != "NPS":
            continue
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
                continue
        temps = DisplayLoggedTemperature.temperatures_from_file(core_file)
        core_temperatures = DisplayLoggedTemperature.average_temperatures(temps)
        date = core_temperatures[0][0].date()
        start_time = datetime.datetime.combine(date, mouse_datum.start[trial_type].time())
        times, temperatures, baseline_temp = create_relative_temperatures(core_temperatures, start_time, relative_temp=True)
        delta_temperature_at = {0: 0}
        average_temperature = None
        for i, time in enumerate(times):
            if time > 0 and not average_temperature:
                average_temperature = np.mean(np.array(temperatures[i:])) + baseline_temp
            if time < .75:
                delta_temperature_at[0.75] = temperatures[i]
            if time < 4:
                delta_temperature_at[4] = temperatures[i]
            if time < 5:
                delta_temperature_at[5] = temperatures[i]
        k = 0
        print("%s:%s:%.3f" % (mouse_number, mouse_datum.vector, average_temperature))
        for k, time in enumerate(times):
            if time > 0:
                break
        f.write(",".join([mouse_number, mouse_datum.genotype, mouse_datum.vector] +
                         ["%.2f" % (s + baseline_temp) for s in delta_temperature_at.values()]))
        f.write(os.linesep)
        print(mouse_number, mouse_datum.genotype, mouse_datum.vector)
        mouse_type = "%s %s" % (mouse_datum.genotype, mouse_datum.vector)
        if mouse_type not in baseline_temps.keys():
            baseline_temps[mouse_type] = []
        baseline_temps[mouse_type].append((average_temperature, delta_temperature_at[0.75], delta_temperature_at[4]))
        line, = ax[1].plot(times, temperatures, c=mouse_datum.get_color(), linestyle=mouse_datum.get_linestyle(), label=mouse_datum.number)
        line_event_handlers[1].add_line(line)
        fig.canvas.mpl_connect("motion_notify_event", line_event_handlers[1].hover)
        times, temperatures, _ = create_relative_temperatures(env_temperatures, start_time, relative_temp=False)
        line, = ax[0].plot(times, temperatures, c=mouse_datum.get_color(), linestyle=mouse_datum.get_linestyle(), label=mouse_datum.number)
        line_event_handlers[0].add_line(line)
        fig.canvas.mpl_connect("motion_notify_event", line_event_handlers[0].hover)
    f.close()
    ax[1].set_xlim([-1, 5])
    ax[1].set_ylim([-20, 5])
    ax[1].axvline(x=0, c='b')
    ax[1].set_ylabel("ΔCore Temp (C)")
    ax[0].set_xlim([-1, 5])
    ax[0].axvline(x=0, c='b')
    ax[0].set_ylabel("Cage Temp (C)")
    ax[1].set_xlabel("Hours from trial start")
    plt.tight_layout()
    plt.show()
    fig, ax = plt.subplots()
    wobble_plot(baseline_temps, ax, lambda x: x[0])
    ax.set_ylabel("Average core temperature during test")
    plt.show()

    fig, ax = plt.subplots()
    wobble_plot(baseline_temps, ax, lambda x: x[2])
    ax.set_ylabel("ΔCore temperature (C) over 4 hours")
    plt.show()
    fig, ax = plt.subplots()
    wobble_plot(baseline_temps, ax, lambda x: x[1])
    ax.set_ylabel("ΔCore temperature (C) over 45 minutes")
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
    plot_condition(mouse_numbers, TrialType.Cold)
    plot_condition(mouse_numbers, TrialType.Hot)
    plot_condition(mouse_numbers, TrialType.Control)