import os
import matplotlib.pyplot as plt
from matplotlib import font_manager
import glob
import scipy.stats
import numpy as np


class Mouse:
    def __init__(self, mouse_number):
        self.mouse_number = mouse_number
        self.trials = []
        self.temperatures, self.hist, self.bin_edges, self.mean, self.bin_centers = None, None, None, None, None

    def calculate_mean_preferences(self):
        self.temperatures = np.concatenate([trial.temperatures for trial in self.trials])
        self.hist, self.bin_edges = np.histogram(self.temperatures, bins=50, range=(0, 50))
        self.mean = np.mean(self.temperatures)
        self.bin_centers = np.zeros(shape=(self.bin_edges.size - 1,))
        for i in range(self.bin_centers.size):
            self.bin_centers[i] = (self.bin_edges[i] + self.bin_edges[i + 1]) / 2

    def plot(self, ax, x, facecolor, edgecolor):
        x1 = np.array([x] * self.bin_centers.size)
        x2 = x1 + (self.hist / (0.3 * np.sum(self.hist)))
        ax.fill_betweenx(self.bin_centers, x1, x2, color=facecolor)
        ax.plot(x2, self.bin_centers, color=edgecolor)
        ax.text(x, self.bin_centers[-1] - x % 2, self.mouse_number)


class Trial:
    def __init__(self, temperatures):
        self.temperatures = temperatures
        self.hist, self.bin_edges = np.histogram(temperatures, bins=50, range=(0, 50))
        self.mean = np.mean(temperatures)
        self.mouse_number = mouse_number
        self.bin_centers = np.zeros(shape=(self.bin_edges.size - 1,))
        for i in range(self.bin_centers.size):
            self.bin_centers[i] = (self.bin_edges[i] + self.bin_edges[i + 1]) / 2

    def plot(self, ax, x, facecolor, edgecolor):
        x1 = np.array([x] * self.bin_centers.size)
        x2 = x1 + (self.hist / 300)
        ax.fill_betweenx(self.bin_centers, x1, x2, color=facecolor)
        ax.plot(x2, self.bin_centers, color=edgecolor)
        ax.text(x, self.bin_centers[-1] - x % 2, self.mouse_number)


def read_keyfile(keyfile_path):
    mice = []
    with open(keyfile_path) as f:
        for line in f:
            mice.append(line.strip().split(":"))
    return mice


def bar_graph(mCh_mice, Cas3_mice, func, ylabel=None):
    mCh_values = []
    for mouse in mCh_mice:
        values = []
        for trial in mouse.trials:
            values.append(func(trial))
        mCh_values.append(sum(values) / len(values))
    average = sum(mCh_values) / len(mCh_values)
    plt.bar([1], [average], color="lightgray")
    plt.scatter([1] * len(mCh_values), mCh_values, facecolor='#b0b0b0', edgecolor='k')

    Cas3_values = []
    for mouse in Cas3_mice:
        values = []
        for trial in mouse.trials:
            values.append(func(trial))
        Cas3_values.append(sum(values) / len(values))
    average = sum(Cas3_values) / len(Cas3_values)
    plt.bar([2], [average], color="#ffc7b3")
    plt.scatter([2] * len(Cas3_values), Cas3_values, facecolor='orangered', edgecolor='k')

    plt.xlim([0.5, 2.5])
    plt.xticks([1, 2], ["mCh", "Cas3"])
    if ylabel:
        plt.ylabel(ylabel)

    _, p = scipy.stats.ttest_ind(mCh_values, Cas3_values)
    plt.title("p=%.3f" % p)
    plt.show()


if __name__ == "__main__":
    font_path = r"C:\Users\Fillan\Documents\Fonts\Myriad Pro Regular.ttf"
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = prop.get_name()
    plt.rcParams['font.size'] = 14
    plt.rcParams['axes.linewidth'] = 2

    keyfile_path = r"R:\Fillan\Parabrachial Ablations\Temperature Bridge\Keyfile.txt"
    mice = read_keyfile(keyfile_path)
    folder = os.path.dirname(keyfile_path)
    time_range = [30, 150]

    Cas3_mice, mCh_mice = [], []
    for mouse_number, vector in mice:
        tsv_files = glob.glob(os.path.join(folder, "*", "%s.tsv" % mouse_number))
        mouse = Mouse(mouse_number)
        for tsv_file in tsv_files:
            temperatures = []
            with open(tsv_file) as f:
                for line in f:
                    time, temperature = line.strip().split("\t")
                    time, temperature = float(time), float(temperature)
                    if time_range[0] < time < time_range[1]:
                        temperatures.append(temperature)
            temperatures = np.array(temperatures)
            mouse.trials.append(Trial(temperatures))
        if vector == "Cas3":
            Cas3_mice.append(mouse)
        if vector == "mCh":
            mCh_mice.append(mouse)

    Cas3_mice = sorted(Cas3_mice, key=lambda x: x.mouse_number)
    mCh_mice = sorted(mCh_mice, key=lambda x: x.mouse_number)
    x_position = 5
    fig, ax = plt.subplots()
    for mouse in Cas3_mice:
        for trial in mouse.trials:
            trial.plot(ax, x_position, edgecolor='orangered', facecolor="#ffc7b3")
            x_position += 1
        x_position += 2
    Cas3_mean = x_position / 2 + 0.5
    x_position += 2
    for mouse in mCh_mice:
        for trial in mouse.trials:
            trial.plot(ax, x_position, edgecolor='#b0b0b0', facecolor="lightgray")
            x_position += 1
        x_position += 2
    mCh_mean = (Cas3_mean + x_position) / 2 + 6.5
    plt.xticks([Cas3_mean, mCh_mean], ["Cas3", "mCh"])
    plt.ylabel("Preferred Temperature (°C)")
    plt.show()

    x_position = 2
    fig, ax = plt.subplots()
    for mouse in Cas3_mice:
        mouse.calculate_mean_preferences()
        mouse.plot(ax, x_position, edgecolor='orangered', facecolor="#ffc7b3")
        x_position += 1
    Cas3_mean = x_position / 2 + 0.5
    x_position += 2
    for mouse in mCh_mice:
        mouse.calculate_mean_preferences()
        mouse.plot(ax, x_position, edgecolor='#b0b0b0', facecolor="lightgray")
        x_position += 1
    mCh_mean = (Cas3_mean + x_position) / 2 + 2.5
    plt.xticks([Cas3_mean, mCh_mean], ["Cas3", "mCh"])
    plt.ylabel("Preferred Temperature (°C)")
    plt.show()

    fig, ax = plt.subplots()
    hist, bin_centers = np.zeros(shape=(50,)), np.zeros(shape=(50,))
    for mouse in Cas3_mice:
        hist += mouse.hist
        bin_centers = mouse.bin_centers
    x1 = np.array([2] * bin_centers.size)
    x2 = x1 + (7 * hist / np.sum(hist))
    ax.fill_betweenx(bin_centers, x1, x2, color="#ffc7b3")
    ax.plot(x2, bin_centers, color='orangered')

    hist, bin_centers = np.zeros(shape=(50,)), np.zeros(shape=(50,))
    for mouse in mCh_mice:
        hist += mouse.hist
        bin_centers = mouse.bin_centers
    x1 = np.array([3] * bin_centers.size)
    x2 = x1 + (7 * hist / np.sum(hist))
    ax.fill_betweenx(bin_centers, x1, x2, color="lightgray")
    ax.plot(x2, bin_centers, color='#b0b0b0')
    plt.xticks([2.5, 3.5], ["Cas3", "mCh"])
    plt.show()

    bar_graph(mCh_mice, Cas3_mice, lambda x: x.mean, "Mean Temperature")
    min_temperature = 25
    bar_graph(mCh_mice, Cas3_mice, lambda x: np.sum(x.temperatures < min_temperature) / 6, "Minutes below %i °C" % min_temperature)
    max_temperature = 45
    for max_temperature in range(35, 45):
        bar_graph(mCh_mice, Cas3_mice, lambda x: np.sum(x.temperatures > max_temperature) / 6, "Minutes above %i °C" % max_temperature)