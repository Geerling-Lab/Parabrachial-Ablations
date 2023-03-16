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
    print(os.path.join(movement_file_directory, "*.tsv"), glob.glob(os.path.join(movement_file_directory, "*.tsv")))
    for file in glob.glob(os.path.join(movement_file_directory, "*.tsv")):
        mouse_number = os.path.basename(file)[:-4]
        print(mouse_number)
        datetimes = []
        movements = []
        first_datetime = None
        last_datetime = None
        lights_on_time = datetime.time(hour=6)
        lights_off_time = datetime.time(hour=18)#MouseLookup(mouse_number, ["Light on", "Light off"])
        assert type(lights_on_time) == datetime.time
        assert type(lights_off_time) == datetime.time
        date1 = datetime.date(1, 1, 1)
        date2 = datetime.date(1, 1, 2)
        lights_off_duration = datetime.datetime.combine(date2, lights_on_time) - \
                              datetime.datetime.combine(date1, lights_off_time)
        with open(file) as f:
            for line in f:
                try:
                    dt, movement = line.strip().split("\t")
                    movement = int(movement)
                except ValueError:
                    continue
                dt_datetime = datetime.datetime.strptime(dt, "%m/%d/%Y %H:%M:%S")
                """
                If we're too close to when the lights turn on or off, there will be a huge amount of movement in the frame
                Therefore, if we're within 90 seconds of the light on or light off time, throw out the frame
                """
                if (datetime.datetime.combine(dt_datetime.date(), lights_off_time) - dt_datetime).seconds < 90:
                    movement = 0
                if (dt_datetime - datetime.datetime.combine(dt_datetime.date(), lights_off_time)).seconds < 90:
                    movement = 0
                if (datetime.datetime.combine(dt_datetime.date(), lights_on_time) - dt_datetime).seconds < 90:
                    movement = 0
                if (dt_datetime - datetime.datetime.combine(dt_datetime.date(), lights_on_time)).seconds < 90:
                    movement = 0
                movements.append(movement)
                if first_datetime is None:
                    first_datetime = dt_datetime
                last_datetime = dt_datetime
                dt_int = time.mktime(dt_datetime.timetuple())
                datetimes.append(dt_int)

        datetimes = np.array(datetimes)
        movements = np.array(movements)
        threshold = 1750  # Chosen by FG and JCG 3/25/21
        percent_epochs_with_movement = np.sum(movements > threshold) / np.size(movements)
        """
        Then, we plot it as a stem plot. We had to convert the x axis to ints from datetimes
        """
        fig, ax = plt.subplots()
        plt.stem(datetimes, movements, use_line_collection=True, markerfmt="None", linefmt="k", basefmt="none")
        plt.title("%s: %.1f%%" % (mouse_number, percent_epochs_with_movement * 100))


        """
        Then, lots of math to get it to display labels every six hours on the X axis
        """
        timedelta = datetime.timedelta(hours=4)
        num_timedeltas = math.ceil((last_datetime - first_datetime) / timedelta)
        first_label_datetime = datetime.datetime(year=first_datetime.year, month=first_datetime.month,
                                                 day=first_datetime.day, hour=first_datetime.hour // 6 * 6)
        label_datetime_strs = []
        label_datetime_floats = []
        for i in range(num_timedeltas + 1):
            dt = first_label_datetime + timedelta * i
            label_datetime_strs.append(dt.strftime("%m/%d\n%H:%M"))
            label_datetime_floats.append(time.mktime(dt.timetuple()))
        plt.xticks(label_datetime_floats, label_datetime_strs)
        #plt.yticks([])
        plt.xlabel("Time")
        plt.ylabel("Movement")
        """
        We make grey areas during lights-off
        """
        plt.xlim(
            time.mktime(first_datetime.timetuple()),
            time.mktime(last_datetime.timetuple()))
        xmin, xmax, ymin, ymax = plt.axis()
        #plt.ylim([0, 9e4])
        i = 0
        while True:
            lights_off = datetime.datetime.combine(first_datetime.date(), lights_off_time) + (i - 1) * datetime.timedelta(days=1)
            width = lights_off_duration.seconds
            if lights_off < first_datetime:
                width -= (first_datetime - lights_off).total_seconds()
                lights_off = first_datetime
            if time.mktime(lights_off.timetuple()) > xmax:
                break
            rect = patches.Rectangle(xy=(time.mktime(lights_off.timetuple()), 0), width=width, height=ymax, ec=None, fc=(0, 0, 1, 0.1))
            ax.add_patch(rect)
            i += 1

        ax.hlines(y=threshold, xmin=xmin, xmax=xmax, color='r', linestyle='-')
        plt.show()
