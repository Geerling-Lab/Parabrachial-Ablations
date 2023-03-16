import numpy as np
import glob
import os
import scipy.integrate, scipy.signal
import matplotlib.pyplot as plt


def moving_average_numpy(t, times, buffer):
    output = np.zeros(np.size(times))
    for i in range(np.size(times)):
        output[i] = moving_average(t, times[i], buffer)
    return output


def moving_average(t, time, buffer):
    time_min = max(0, time - buffer)
    time_max = min(np.max(t[:, 0]), time + buffer)
    return np.average(t[np.logical_and(time_min < t[:, 0], t[:, 0] < time_max), 1])


def to_string(n):
    s = ""
    for i in range(np.size(n)):
        s += (str(n[i]) + ",")
    return s


if __name__ == "__main__":
    mCh_mice = ["02547", "02548", "02549", "02550", "02551", "02905", "02909", "02910", "02910", "03433", "03434"]
    directory = r"/home/fillan/rdss_drive/Fillan/Parabrachial Arousal/MVT"
    file_paths = glob.glob(os.path.join(directory, "* Temperature.csv"))
    with open(os.path.join(directory, "AverageTemperature.csv"), "w+") as f:
        for file_path in file_paths:
            animal_name = os.path.basename(file_path).split(" ")[0]
            t = np.genfromtxt(file_path, delimiter=",")[1:, 2:]
            x_new = np.arange(0, 8000, 100)
            y_new = moving_average_numpy(t, x_new, 500)
            if animal_name in mCh_mice:
                color = "b-"
            else:
                color = "r-"
            plt.plot(x_new, y_new, color)
            average_temperature = scipy.integrate.trapz(t[:, 1], t[:, 0]) / (np.max(t[:, 0] - np.min(t[:, 0])))
            f.write("%s,%2f,%s\n" % (animal_name, average_temperature, to_string(y_new)))
    plt.show()