import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import datetime
import glob
import os
import time
import numpy as np
from matplotlib import font_manager
#from MouseLookup import MouseLookup


def read_tsv_file(file_path, approximate_lights_on, approximate_lights_off):
    datetimes = []
    movements = []
    lights_off_times = []
    lights_on_times = []
    first_line = True
    with open(file_path) as f:
        for line in f:
            try:
                dt, movement = line.strip().split("\t")
                movement = int(movement)
            except ValueError:
                continue
            if first_line:
                first_line = False
                continue
            dt_datetime = datetime.datetime.strptime(dt, "%m/%d/%Y %H:%M:%S")
            """
            If we're too close to when the lights turn on or off, there will be a huge amount of movement in the frame
            """
            if movement > 2e4 and abs((dt_datetime - datetime.datetime.combine(dt_datetime.date(),
                                                                           approximate_lights_on)).total_seconds()) < 15 * 60:
                print("Lights on time: %s" % dt_datetime)
                movement = 0
                lights_on_times.append(dt_datetime.time())
            if movement > 2e4 and abs((dt_datetime - datetime.datetime.combine(dt_datetime.date(),
                                                                           approximate_lights_off)).total_seconds()) < 15 * 60:
                print("Lights off time: %s" % dt_datetime)
                movement = 0
                lights_off_times.append(dt_datetime.time())
            movements.append(movement)
            datetimes.append(dt_datetime)

    if len(lights_on_times) == 0:
        lights_on_time = approximate_lights_on
    else:
        average_seconds = int(sum([int(x.hour * 3600 + x.minute * 60 + x.second) for x in lights_on_times]) / len(lights_on_times))
        hour = average_seconds // 3600
        minute = (average_seconds % 3600) // 60
        second = average_seconds % 60
        lights_on_time = datetime.time(hour=hour, minute=minute, second=second)

    if len(lights_off_times) == 0:
        lights_off_time = approximate_lights_off
    else:
        average_seconds = int(sum([int(x.hour * 3600 + x.minute * 60 + x.second) for x in lights_off_times]) / len(lights_off_times))
        hour = average_seconds // 3600
        minute = (average_seconds % 3600) // 60
        second = average_seconds % 60
        lights_off_time = datetime.time(hour=hour, minute=minute, second=second)

    return datetimes, movements, lights_on_time, lights_off_time


if __name__ == "__main__":
    """
    Generates running-wheel plots from Movement Files
    Movement Files are produced by the VideoScoring.py program
    Displays the amount of movement (pixels moved) as a function of time
    """
    font_path = r"C:\Users\richi\PycharmProjects\Fonts\Myriad Pro Regular.ttf"
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = prop.get_name()
    plt.rcParams['font.size'] = 18
    plt.rcParams['axes.linewidth'] = 2
    """
    First, we load in the data from the movement file
    """

    movement_file_directory = r"C:\Users\richi\PycharmProjects\VideoScoring\Movement Files"
    for file_path in glob.glob(os.path.join(movement_file_directory, "*.tsv")):
        mouse_number = os.path.basename(file_path)[:-4]
        print(mouse_number)
        if "3453" in mouse_number or "3558" in mouse_number:
            continue
        if "3559" not in mouse_number:
            continue
        approximate_lights_on = datetime.time(hour=7)
        approximate_lights_off = datetime.time(hour=19)
        datetimes, movements, lights_on, lights_off = read_tsv_file(file_path, approximate_lights_on, approximate_lights_off)
        print(lights_on, lights_off)

        date1 = datetime.date(1, 1, 1)
        date2 = datetime.date(1, 1, 2)
        lights_off_duration = datetime.datetime.combine(date2, lights_on) - \
                              datetime.datetime.combine(date1, lights_off)

        datetimes = np.array(datetimes)
        movements = np.array(movements)
        threshold = 1750  # Chosen by FG and JCG 3/25/21
        percent_epochs_with_movement = np.sum(movements > threshold) / np.size(movements)
        """
        Then, we plot it as a stem plot. We had to convert the x axis to ints from datetimes
        """
        fig, axs = plt.subplots(nrows=2, sharex=True, gridspec_kw={'height_ratios': [1, 0.08]})
        plt.subplots_adjust(wspace=0, hspace=0)
        axs[0].stem([time.mktime(x.timetuple()) for x in datetimes], movements, use_line_collection=True, markerfmt="None", linefmt="k", basefmt="none")
        axs[0].set_title("%s: %.1f%%" % (mouse_number, percent_epochs_with_movement * 100))


        """
        Then, lots of math to get it to display labels every six hours on the X axis
        """
        first_datetime = datetime.datetime.combine(datetimes[0].date(), lights_on)
        last_datetime = datetime.datetime.combine(datetimes[-1].date(), lights_on)
        timedelta = datetime.timedelta(hours=6)
        num_timedeltas = math.ceil((last_datetime - first_datetime) / timedelta)
        first_label_datetime = datetime.datetime(year=first_datetime.year, month=first_datetime.month,
                                                 day=first_datetime.day, hour=first_datetime.hour // 6 * 6)
        label_datetime_strs = []
        label_datetime_floats = []
        for i in range(num_timedeltas + 1):
            dt = first_label_datetime + timedelta * i
            label_datetime_strs.append(dt.strftime("%m/%d\n%H:%M"))
            label_datetime_floats.append(time.mktime(dt.timetuple()))
        axs[0].set_xticks(label_datetime_floats, label_datetime_strs)
        axs[0].set_yticks([])
        axs[1].set_xlabel("Time")
        axs[0].set_ylabel("Movement")
        """
        We make grey areas during lights-off
        """
        axs[0].set_xlim(
            time.mktime(first_datetime.timetuple()),
            time.mktime(last_datetime.timetuple()))
        ymin, ymax = axs[0].get_ylim()
        xmin, xmax = axs[0].get_xlim()
        axs[0].set_ylim([0, ymax])
        i = 0
        rect = patches.Rectangle(xy=(xmin, 0), width=(xmax-xmin), height=1, facecolor=(255 / 255, 251 / 255, 36 / 255))
        axs[1].add_patch(rect)
        while True:
            try:
                this_lights_off = datetime.datetime.combine(datetime.date(year=first_datetime.year, month=first_datetime.month,
                                           day=first_datetime.day + i - 1), lights_off)
            except:
                this_lights_off = datetime.datetime.combine(datetime.date(year=first_datetime.year, month=first_datetime.month + 1,
                                                                     day=(first_datetime.day + i - 1) % 30), lights_off)
            width = lights_off_duration.seconds
            if this_lights_off < first_datetime:
                width -= (first_datetime - this_lights_off).total_seconds()
                this_lights_off = first_datetime
            if time.mktime(this_lights_off.timetuple()) > xmax:
                break
            rect = patches.Rectangle(xy=(time.mktime(this_lights_off.timetuple()), 0), width=width, height=1, facecolor=(30 / 255, 11 / 255, 179 / 255))
            axs[1].add_patch(rect)
            i += 1
        axs[1].set_yticks([])
        axs[0].hlines(y=threshold, xmin=xmin, xmax=xmax, color='r', linestyle='-')
        axs[0].hlines(y=1e5, xmin=xmin, xmax=xmax, color='r', linestyle='-')
        plt.show()
