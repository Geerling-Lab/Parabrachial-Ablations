import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines
import glob
import os
import sys
from enum import IntEnum
from datetime import datetime, timedelta
import operator
import pandas as pd
import scipy.interpolate


class SleepState(IntEnum):
    Wake = 1
    NREM = 2
    REM = 3
    Undefined = 4


class ReadPowers:
    """
    This class reads in a frequency file exported by Sirenia Sleep Pro
    The frequency file is split up into the frequencies, hour by hour, as well as by sleep state
    So, for example, the frequency distribution of NREM sleep from 8:00-9:00 is reported
    To average these together, we weigh the frequency distribution by how much of each state (i.e. NREM) there was
    Therefore, this requires an instance of ReadScores()
    """
    def __init__(self, file_name, rs):
        self.file_name = file_name
        self.number = os.path.basename(os.path.splitext(file_name)[0])
        self.rs = rs
        self.minimum_time = datetime(2040, 1, 1, 1, 1, 1, 1)
        self.maximum_time = datetime(1900, 1, 1, 1, 1, 1, 1)
        df = pd.read_csv(self.file_name, header=5, engine='python', delimiter="\t", error_bad_lines=False,
                         warn_bad_lines=False)
        self.powers = np.zeros(shape=(5, df.shape[0]), dtype=np.float)
        one_hour = timedelta(hours=1)
        """
        self.powers[0,:] is the frequencies
        self.powers[1,:] is the power of Wake at each frequency
        self.powers[2,:] is the power of NREM at each frequency
        self.powers[3,:] is the power of REM at each frequency
        as defined by int(SleepState)
        """
        self.powers[0, :] = df["Time from Start"]
        for col in df.columns:
            if "Other" in col:
                ss = SleepState.Undefined
            elif "Wake" in col:
                ss = SleepState.Wake
            elif "Non REM" in col:
                ss = SleepState.NREM
            elif "REM" in col:
                ss = SleepState.REM
            else:
                continue
            try:  # If any of the entries are not formatted as floats or >1e6, then throw out the entire column
                if np.sum(df[col] > 1e6) > 0:
                    continue
            except TypeError as te:
                continue
            dt_string = " ".join(col.split(" ")[-2:])
            dt = datetime.strptime(dt_string, "%m/%d/%y %H:%M:%S")
            if dt < self.minimum_time:
                self.minimum_time = dt
            if dt > self.maximum_time:
                self.maximum_time = dt
            if self.rs is not None:
                self.powers[int(ss), :] += \
                    df[col] * rs.time_in_state(start_date_time=dt, end_date_time=dt+one_hour, sleep_state=ss)
            else:
                self.powers[int(ss), :] += df[col]
        for i in range(1, 5):
            if self.rs is not None:
                self.powers[i, :] /= rs.totals[i]
            else:
                pass

    def graph(self, ax):
        ax.plot(self.powers[0, :], self.powers[1, :], "r-")  # Wake
        ax.plot(self.powers[0, :], self.powers[2, :], "g-")  # NREM
        ax.plot(self.powers[0, :], self.powers[3, :], "b-")  # REM
        ax.plot(self.powers[0, :], self.powers[4, :], "k-")
        ax.set_xticklabels([])
        ax.set_xlim([0.5, 20])
        ax.set_ylabel(self.number)


class ReadScores():
    """
    This class reads in a score file exported by Sirenia Sleep Pro
    It is only designed to be accessed via the "time_in_state" function
    """
    def __init__(self, file_name):
        self.scores = []
        self.totals = [0, 0, 0, 0, 0]
        self.minimum_time = datetime(2040, 1, 1, 1, 1, 1, 1)
        self.maximum_time = datetime(1900, 1, 1, 1, 1, 1, 1)
        with open(file_name) as f:
            for line in f:
                words = line.split("\t")
                if len(words) >= 5:
                    try:
                        date_time_string = "%s %s" % (words[0].strip(), words[1].strip())
                        date_time = datetime.strptime(date_time_string, "%m/%d/%Y %H:%M:%S")
                        if date_time > self.maximum_time:
                            self.maximum_time = date_time
                        if date_time < self.minimum_time:
                            self.minimum_time = date_time
                        sleep_state_number = int(float(words[4].strip()))
                        if sleep_state_number == 1:
                            sleep_state = SleepState.Wake
                        elif sleep_state_number == 2:
                            sleep_state = SleepState.NREM
                        elif sleep_state_number == 3:
                            sleep_state = SleepState.REM
                        else:
                            sleep_state = SleepState.Undefined
                        self.totals[int(sleep_state)] += 1
                        self.scores.append((date_time, sleep_state))
                    except ValueError as ve:
                        pass
        self.scores.sort(key=operator.itemgetter(0))

    def time_in_state(self, start_date_time, end_date_time, sleep_state):
        in_state = 0
        for date_time, ss in self.scores:
            if start_date_time < date_time < end_date_time and sleep_state == ss:
                in_state += 1
        return in_state


class AveragePowers:
    def __init__(self):
        self.powers = np.zeros((5, 300))
        self.powers[0, :] = np.arange(0, 30, 0.1)

    def add_mouse(self, powers):
        for i in range(1, 5):
            x = np.nan_to_num(powers[0, :])
            if np.max(powers[i, :]) == 0:
                y = powers[i, :]
            else:
                y = np.nan_to_num(powers[i, :] / np.max(powers[i, :]))
            f = scipy.interpolate.interp1d(x, y, kind="cubic")
            self.powers[i, :] += f(self.powers[0, :])

    def normalize(self, num_mice):
        for i in range(1, 5):
            self.powers[i, :] /= num_mice

    def save(self, file_path):
        np.savetxt(file_path, self.powers, delimiter=",")

    def __getitem__(self, item):
        return self.powers[item]


if __name__ == "__main__":
    """
    rs = ReadScores(r"R:\\Shantelle\\TSV Export\\02549.tsv")
    rp = ReadPowers(r"R:/Shantelle/Frequency Export/2549.txt", rs)
    fig, axes = plt.subplots(2, 2)
    rp.graph(axes[0, 0])
    plt.show()
    """
    Master = {}
    colors = ['b', 'r', 'g', 'y']
    power_folder = r"R:\Shantelle\Frequency Export"
    score_folder = r"R:\Shantelle\TSV Export"
    with open(os.path.join(power_folder, "Master.txt")) as f:
        for line in f:
            number, histology = line.strip().split(" ")
            if histology in Master:
                Master[histology].append(number)
            else:
                Master[histology] = [number]
    max_length = 0
    for key, item in Master.items():
        if len(item) > max_length:
            max_length = len(item)
    fig, axes = plt.subplots(ncols=len(Master.items()), nrows=max_length+1)
    count = 0
    first_frequency = None
    for genotype, mouse_list in Master.items():
        average_powers = AveragePowers()
        for i, mouse_number in enumerate(mouse_list):
            power_file = os.path.join(power_folder, "%04d.txt" % int(mouse_number))
            score_file = os.path.join(score_folder, "%04d.tsv" % int(mouse_number))
            if not os.path.exists(power_file):
                continue
            if not os.path.exists(score_file):
                rs = None
            else:
                rs = ReadScores(score_file)
            rp = ReadPowers(power_file, rs)
            """
            if rs.maximum_time < rp.minimum_time or rs.minimum_time > rp.maximum_time:
                print("Times for %s don't line up" % mouse_number)
                print("Score file runs from %s to %s" % (rs.minimum_time, rs.maximum_time))
                print("Power file runs from %s to %s" % (rp.minimum_time, rp.maximum_time))
            """
            rp.graph(axes[i+1, count])
            average_powers.add_mouse(rp.powers)
        average_powers.normalize(len(mouse_list))
        average_powers.save("%s.csv" % genotype)
        axes[0, count].plot(average_powers[0, :], average_powers[1, :], "r-")
        axes[0, count].plot(average_powers[0, :], average_powers[2, :], "g-")
        axes[0, count].plot(average_powers[0, :], average_powers[3, :], "b-")
        axes[0, count].plot(average_powers[0, :], average_powers[4, :], "k-")
        axes[0, count].set_xticklabels([])
        axes[0, count].set_xlim([0.5, 20])
        axes[0, count].set_ylabel("Average")
        count += 1
    axes[6, 0].set_xlim([0.5, 20])
    red_line = matplotlib.lines.Line2D([], [], color="red", label="Wake")
    blue_line = matplotlib.lines.Line2D([], [], color="green", label="NREM")
    green_line = matplotlib.lines.Line2D([], [], color="blue", label="REM")
    axes[6, 0].legend(handles=[red_line, blue_line, green_line], loc="upper right")
    plt.show()