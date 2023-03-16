import os.path
import matplotlib.pyplot as plt
import numpy as np
import datetime
import DisplayLoggedTemperature
import statistics
from matplotlib import font_manager
import scipy.stats

PLUNGE_TEMPERATURES = range(30, 20, -1)


def find_cold_start(movement_file):
    """
    :param movement_file:
    :return:
    This is a helper function that finds breaks in a movement file
    I can determine when mice are placed in the fridge by the break in the movement file (as created by vLMA.py)
    However, occasionally there is more than one break. This function displays all their breaks to allow
    you manual selection of which break corresponds to putting the mouse in the fridge
    """
    interval_begin = None
    last_time = None
    with open(movement_file) as f:
        breaks = []
        for line in f:
            line = line.strip()
            words = line.split(",")
            if len(words) != 2:
                print(line)
            time, _ = words
            time = datetime.datetime.strptime(time, "%m/%d/%y %H:%M:%S")
            if not interval_begin:
                interval_begin = time
            else:
                if last_time:
                    time_delta = time - last_time
                    if time_delta.total_seconds() > 20:
                        breaks.append((interval_begin, last_time))
                        interval_begin = time
                last_time = time
    breaks.append((interval_begin, last_time))
    for i in range(len(breaks) - 1):
        first_interval_length = breaks[i][1] - breaks[i][0]
        second_interval_length = breaks[i+1][1] - breaks[i+1][0]
        break_length = breaks[i+1][0] - breaks[i][1]
        break_time = breaks[i][1]
        print("%s\t\t%s\t\t%s\t\t%s" % (first_interval_length, second_interval_length, break_length, break_time))


def get_color(Cas3, fed):
    if Cas3:
        if fed:
            return "#FFB19F"
        else:
            return "orangered"
    else:
        if fed:
            return "lightgray"
        else:
            return "black"


class Condition:
    def __init__(self, Cas3, food):
        self.Cas3 = Cas3
        self.food = food
        self.cold_times = {}
        self.movements = []
        self.bats = []
        self.tails = []
        self.mean_temperature = []

class Mouse:
    def __init__(self, mouse_number, Cas3):
        self.mouse_number = mouse_number
        self.Cas3 = Cas3
        self.axs = None
        self.trials = []

    def graph_individual(self):
        for trial in self.trials:
            color = get_color(self.Cas3, trial.fed)
            self.axs[0].plot(trial.core_times, trial.core_temperatures, color=color)
            self.axs[1].plot(trial.movement_times, trial.movements, color=color)
            self.axs[2].plot(trial.bat_times, trial.bat_temperatures, color=color)
            self.axs[3].plot(trial.tail_times, trial.tail_temperatures, color=color)
            self.axs[4].plot(trial.environment_times, trial.environment_temperatures, color=color)

    @staticmethod
    def get_average_between_times(times, values, start_time, end_time, carry_last_point=False):
        # start_time and end_time should be expressed as hours (float), not datetimes
        start_index = None
        end_index = None
        if len(values) == 0:
            return None
        assert len(times) == len(values)
        for i, time in enumerate(times):
            if time > start_time and start_index is None:
                start_index = i
            if time > end_time:
                end_index = i
                break
        if start_index is None and end_index is None:
            if carry_last_point:
                return values[-1]
            else:
                return None
        included_values = values[start_index: end_index]
        if len(included_values) > 0:
            return sum(included_values) / len(included_values)
        else:
            assert start_index == end_index
            start_index -= 1
            middle_time = (start_time + end_time) / 2
            return values[start_index] + (values[end_index] - values[start_index])/(times[end_index] - times[start_index]) * (middle_time - times[start_index])

    @staticmethod
    def get_value_at_time(times, values, time):
        if time < times[0]:
            return values[0]
        for i in range(len(times)):
            if times[i] > time:
                return values[i-1] + (values[i] - values[i-1]) / (times[i] - times[i-1]) * (time - times[i-1])
        else:
            return values[-1]


class Trial:
    def __init__(self, fed, folder, fridge_time, mouse_number):
        self.fed = fed
        self.mouse_number = mouse_number
        self.fridge_time = fridge_time
        self.movement_file = os.path.join(folder, "%s Movement.txt" % mouse_number)
        if not os.path.exists(self.movement_file):
            raise AssertionError("%s does not exist" % self.movement_file)
        self.core_file = os.path.join(folder, "%s Core.csv" % mouse_number)
        if not os.path.exists(self.core_file):
            raise AssertionError("%s does not exist" % self.core_file)
        self.flir_file = os.path.join(folder, mouse_number, "Output.csv")
        if not os.path.exists(self.flir_file):
            raise AssertionError("%s does not exist" % self.flir_file)
        self.rump_file = os.path.join(folder, mouse_number, "Rump.csv")
        if not os.path.exists(self.rump_file):
            pass
            #raise AssertionError("%s does not exist" % self.rump_file)
        self.environmental_file = r"R:\Fillan\Parabrachial Ablations\Temperature in Fridge\Fridge Temperature.csv"  # Downloaded from SensorPush app

        self.movements = []
        self.movement_times = []
        self.core_temperatures = []
        self.core_times = []
        self.tail_temperatures = []
        self.tail_times = []
        self.bat_temperatures = []
        self.bat_times = []
        self.environment_temperatures = []
        self.environment_times = []
        self.rump_temperatures = []
        self.rump_times = []
        self.core_file_begin_time = None
        self.core_file_end_time = None
        self.time_till_cold = {}

        self.read_movement_file()
        self.read_core_file()
        self.read_environmental_file()
        self.read_rump_file()
        self.read_flir_file()

    def read_movement_file(self):
        largest_movement = 0
        with open(self.movement_file) as f:
            for line in f:
                time, movement = line.strip().split(",")
                self.movement_times.append((datetime.datetime.strptime(time, "%m/%d/%y %H:%M:%S") - self.fridge_time).total_seconds() / 3600)
                movement = float(movement)
                if movement > 4 * largest_movement and len(self.movements) > 100:
                    print("Spurious movement: %s at %s" % (self.movement_file, time))
                elif movement > largest_movement:
                    largest_movement = movement
                self.movements.append(movement)

    def read_core_file(self):
        temps = DisplayLoggedTemperature.temperatures_from_file(self.core_file)
        core_temperatures = DisplayLoggedTemperature.average_temperatures(temps)
        for time, temperature in core_temperatures:
            self.core_times.append((time - self.fridge_time).total_seconds() / 3600)
            self.core_temperatures.append(temperature)
            if len(self.time_till_cold) < len(PLUNGE_TEMPERATURES):
                next_plunge_temperature = PLUNGE_TEMPERATURES[(len(self.time_till_cold))]
                if temperature < next_plunge_temperature:
                    self.time_till_cold[next_plunge_temperature] = min((time - self.fridge_time).total_seconds() / 3600, 24)
        for i in range(len(self.time_till_cold), len(PLUNGE_TEMPERATURES)):
            self.time_till_cold[PLUNGE_TEMPERATURES[i]] = 24
        self.core_file_begin_time = self.fridge_time + datetime.timedelta(hours=self.core_times[0])
        self.core_file_end_time = self.fridge_time + datetime.timedelta(hours=self.core_times[-1])

    def read_environmental_file(self):
        with open(self.environmental_file) as f:
            reading_header = True
            for line in f:
                if reading_header:
                    reading_header = False
                    continue
                time, temperature, humidity = line.split(",")
                time = datetime.datetime.strptime(time[1:-1], "%Y-%m-%d %H:%M")
                if not self.core_file_begin_time < time < self.core_file_end_time:
                    continue
                self.environment_times.append((time - self.fridge_time).total_seconds() / 3600)
                if time < self.fridge_time:
                    self.environment_temperatures.append(22)
                else:
                    self.environment_temperatures.append(float(temperature[1:-1]))

    def read_flir_file(self):
        continuously_excluded = 0
        with open(self.flir_file) as f:
            for line in f:
                words = line.strip().split(",")
                time = " ".join(words[0].split(" ")[1:])
                time = (datetime.datetime.strptime(time, "%m %d %Y %H %M %S") - self.fridge_time).total_seconds() / 3600
                bat_temperature = float(words[1])
                if len(self.bat_temperatures) > 10:
                    historical = np.array(self.bat_temperatures[-10:])
                    Z_score = abs((bat_temperature - np.mean(historical)) / np.std(historical))
                    if (Z_score < 6 or continuously_excluded > 3) and bat_temperature < 45:
                        continuously_excluded = 0
                        self.bat_times.append(time)
                        self.bat_temperatures.append(bat_temperature - Mouse.get_value_at_time(self.rump_times, self.rump_temperatures, time))
                    else:
                        continuously_excluded += 1
                else:
                    self.bat_times.append(time)
                    self.bat_temperatures.append(bat_temperature - Mouse.get_value_at_time(self.rump_times, self.rump_temperatures, time))
                if len(words) > 2:
                    if len(words[2]) > 1:
                        tail_temperature = float(words[2])
                        self.tail_times.append(time)
                        self.tail_temperatures.append(tail_temperature)
        self.bat_temperatures = Trial.median_filter(self.bat_temperatures, kernel=5)

    def read_rump_file(self):
        with open(self.rump_file) as f:
            for line in f:
                words = line.strip().split(",")
                if len(words) == 2:
                    try:
                        time = " ".join(words[0].split(" ")[1:])
                        time = (datetime.datetime.strptime(time,
                                                           "%m %d %Y %H %M %S") - self.fridge_time).total_seconds() / 3600
                        rump_temperature = float(words[1])
                        self.rump_temperatures.append(rump_temperature)
                        self.rump_times.append(time)
                    except ValueError as ve:
                        print("%s has an invalid line" % self.rump_file)

    @staticmethod
    def median_filter(list, kernel):
        step = int((kernel - 1) / 2)
        new_list = []
        for i in range(len(list)):
            new_list.append(statistics.median(list[max(0, i - step): min(len(list), i + step + 1)]))
        return new_list


def error_bar(ax, x, y, yerr, color, label=None):
    ax.plot(x, y, c=color, label=label)
    lower = [i - j for i, j in zip(y, yerr)]
    upper = [i + j for i, j in zip(y, yerr)]
    ax.fill_between(x, lower, upper, color=color, alpha=0.1, edgecolor='w')
    #ax.plot(x, lower, c=color, alpha=0.3)
    #ax.plot(x, upper, c=color, alpha=0.3)


def graph_average(mice, ax, ylabel, func_time, func_value, time_delta):
    times = list(np.arange(-24, 24, time_delta))
    values = {}
    """
    values holds all the values (temperatures, movements, ....) for all the mice
    It's organized so that values[(Cas3, fed)] = [[values at first timepoint...], [values at second timepoint...]...]
    """
    for mouse in mice:
        for trial in mouse.trials:
            if (mouse.Cas3, trial.fed) not in values:
                values[(mouse.Cas3, trial.fed)] = [[] for _ in range(len(times))]
            for i, time in enumerate(times):
                start_time = time - time_delta / 2
                end_time = time + time_delta / 2
                average_value = Mouse.get_average_between_times(func_time(trial), func_value(trial), start_time, end_time, carry_last_point=(ylabel != "Environmental\nTemperature"))
                if average_value is not None:
                    values[(mouse.Cas3, trial.fed)][i].append(average_value)
    for (Cas3, fed), value in values.items():
        y, y_err = [], []
        for value_at_timepoint in value:
            value_at_timepoint = np.array(value_at_timepoint)
            y.append(np.mean(value_at_timepoint))
            y_err.append(np.std(value_at_timepoint))
        error_bar(ax, times, y, y_err, color=get_color(Cas3, fed))
        ax.set_ylabel(ylabel)


def plot_condition(conditions, func, ylabel):
    labels = []
    num_statistically_significant = 0
    fig, axs = plt.subplots(nrows=2, gridspec_kw={"height_ratios": [1, 6]})

    conditions_list = list(conditions.values())
    print(ylabel)
    for row_index, condition in enumerate(conditions_list):
        Cas3, food = condition.Cas3, condition.food
        values = [x for x in func(condition) if x is not None]
        xs = [row_index] * len(values)
        plt.scatter(xs, values, c=get_color(Cas3, food))
        average = sum(values) / len(values)
        std = np.std(np.array(values))
        print("%s Cas3=%s, food=%s, %s" % (ylabel, Cas3, food, values))
        plt.plot([row_index - .2, row_index + .2], [average, average], "k-")
        labels.append("%s\n%s" % ("Cas3" if Cas3 else "mCh", "Food" if food else "No Food"))
        for j in range(row_index + 1, len(conditions_list)):
            new_values = func(conditions_list[j])
            _, p = scipy.stats.ttest_ind(values, new_values, equal_var=False)
            if p < 0.1:
                axs[0].plot([row_index, row_index],
                           [num_statistically_significant, num_statistically_significant + .5],
                           "k-")
                axs[0].plot([j, j], [num_statistically_significant, num_statistically_significant + .5],
                           "k-")
                axs[0].plot([row_index, j],
                           [num_statistically_significant + .5, num_statistically_significant + .5], "k-")
                axs[0].text(x=(row_index + j) * 0.5, y=num_statistically_significant + .5, ha="center", va="bottom",
                           s="p=%.4f" % p, fontsize=12)
                num_statistically_significant += 2
    axs[0].set_ylim([0, num_statistically_significant + 1])
    axs[1].spines['top'].set_visible(False)
    axs[0].set_xlim(axs[1].get_xlim())
    axs[0].spines['bottom'].set_visible(False)
    axs[0].set_xticks([])
    axs[0].set_yticks([])
    axs[1].set_ylabel(ylabel)
    axs[1].set_xticks(range(len(labels)))
    axs[1].set_xticklabels(labels)
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

    #find_cold_start(r"R:\Fillan\Parabrachial Ablations\Temperature in Fridge\6523 No Food\6523 Movement.txt")
    keyfile = r"R:\Fillan\Parabrachial Ablations\Temperature in Fridge\Keyfile.txt"
    active_mice = []
    all_mice = []
    reading_mice = True

    with open(keyfile) as f:
        for line in f:
            if line[0] == ">":
                reading_mice = False
                fed, folder, fridge_time = line[1:].strip().split(",")
                folder = folder.strip('"')
                fed = (fed == "Fed")
                fridge_time = datetime.datetime.strptime(fridge_time, "%Y-%m-%d %H:%M:%S")
                for mouse in active_mice:
                    mouse.trials.append(Trial(fed=fed, folder=folder, fridge_time=fridge_time, mouse_number=mouse.mouse_number))
            else:
                if not reading_mice:
                    all_mice.extend(active_mice)
                    active_mice = []
                    reading_mice = True
                mouse_number, genotype = line.strip().split(":")
                genotype = (genotype == "Cas3")
                active_mice.append(Mouse(mouse_number, genotype))
        all_mice.extend(active_mice)

    for mouse in all_mice:
        fig, axs = plt.subplots(nrows=5, sharex=True)
        for trial in mouse.trials:
            color = get_color(mouse.Cas3, trial.fed)
            axs[0].plot(trial.core_times, trial.core_temperatures, color=color)
            axs[0].set_ylim([15, 42])
            axs[1].plot(trial.movement_times, trial.movements, color=color)
            axs[1].set_ylim([0, 120000])
            axs[2].plot(trial.bat_times, trial.bat_temperatures, color=color)
            axs[2].set_ylim([5, 30])
            axs[3].plot(trial.tail_times, trial.tail_temperatures, color=color)
            axs[3].set_ylim([5, 30])
            axs[4].plot(trial.environment_times, trial.environment_temperatures, color=color)
            axs[4].set_ylim([3, 25])
            axs[0].set_title(mouse.mouse_number)
            axs[0].set_ylabel("Core Temperature")
            axs[1].set_ylabel("Movement")
            axs[2].set_ylabel("Intrascapular - Rump")
            axs[3].set_ylabel("Tail")
            axs[4].set_ylabel("Environment")
            axs[4].set_xlim([-1, 22])
        plt.show()

    fig, axs = plt.subplots(nrows=5, sharex=True)
    for mouse in all_mice:
        mouse.axs = axs
        for trial in mouse.trials:
            if not trial.fed:
                color = get_color(mouse.Cas3, trial.fed)
                time_delta = 0.1
                times = list(np.arange(-24, 24, time_delta))
                smoothed_times = []
                smoothed_values = []
                for i, time in enumerate(times):
                    start_time = time - time_delta / 2
                    end_time = time + time_delta / 2
                    smoothed_times.append(time)
                    smoothed_values.append(Mouse.get_average_between_times(trial.core_times, trial.core_temperatures, start_time,
                                                                    end_time, carry_last_point=False))
                axs[0].plot(smoothed_times, smoothed_values, color=color, alpha=0.2, linewidth=1.5)
    graph_average(all_mice, axs[0], "Core\nTemperature", lambda x: x.core_times, lambda x: x.core_temperatures, time_delta=.25)
    graph_average(all_mice, axs[1], "Movement", lambda x: x.movement_times, lambda x: x.movements, time_delta=1)
    graph_average(all_mice, axs[2], "Intrascapular - Rump\nTemperature", lambda x: x.bat_times, lambda x: x.bat_temperatures, time_delta=0.25)
    graph_average(all_mice, axs[3], "Tail\nTemperature", lambda x: x.tail_times, lambda x: x.tail_temperatures, time_delta=.5)
    graph_average(all_mice, axs[4], "Environmental\nTemperature", lambda x: x.environment_times, lambda x: x.environment_temperatures, time_delta=.25)

    axs[0].set_ylim([15, 42])
    axs[1].set_ylim([0, 12000])
    axs[3].set_ylim([5, 30])
    axs[4].set_ylim([3, 25])
    axs[4].set_xlim([-1, 22])

    axs[0].plot([0, 0], [10, 10], color=get_color(Cas3=True, fed=True), label="Cas3, Food")
    axs[0].plot([0, 0], [10, 10], color=get_color(Cas3=True, fed=False), label="Cas3, No Food")
    axs[0].plot([0, 0], [10, 10], color=get_color(Cas3=False, fed=True), label="mCh, Food")
    axs[0].plot([0, 0], [10, 10], color=get_color(Cas3=False, fed=False), label="mCh, No Food")
    axs[0].legend()
    axs[1].set_ylabel("Movement")
    axs[3].set_ylabel("Tail\nTemperature")
    axs[4].set_ylabel("Environmental\nTemperature")
    axs[4].set_xlabel("Hours")
    plt.show()

    conditions = {}
    for Cas3 in [False, True]:
        for fed in [True, False]:
            conditions[(Cas3, fed)] = Condition(Cas3, fed)
    for mouse in all_mice:
        for trial in mouse.trials:
            condition = conditions[(mouse.Cas3, trial.fed)]
            condition.movements.append(Mouse.get_average_between_times(trial.movement_times, trial.movements, 0, 3, carry_last_point=False))
            condition.bats.append(Mouse.get_average_between_times(trial.bat_times, trial.bat_temperatures, 0, 1, carry_last_point=False))
            condition.tails.append(Mouse.get_average_between_times(trial.tail_times, trial.tail_temperatures, 0, 3, carry_last_point=False))
            condition.mean_temperature.append(Mouse.get_average_between_times(trial.core_times, trial.core_temperatures, 0, 24, carry_last_point=False))
            for plunge_temperature, time_till_cold in trial.time_till_cold.items():
                if plunge_temperature not in condition.cold_times:
                    condition.cold_times[plunge_temperature] = []
                condition.cold_times[plunge_temperature].append(time_till_cold)
    for plunge_temperature in PLUNGE_TEMPERATURES:
        plot_condition(conditions, lambda x: x.cold_times[plunge_temperature], "Hours until core temperature = %i C" % plunge_temperature)
    plot_condition(conditions, lambda x: x.mean_temperature, "Mean temperature (Cold exposure)")
    plot_condition(conditions, lambda x: x.movements, "Movements (Hours 0-3)")
    plot_condition(conditions, lambda x: x.bats, "Average BAT (Hours 0-3)")
    plot_condition(conditions, lambda x: x.tails, "Average Tail (Hours 0-3)")
