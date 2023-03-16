import argparse
import imutils
import time
from PIL import Image
import cv2
import numpy as np
import datetime
import ctypes
import matplotlib.pyplot as plt
import os
import pandas as pd


class VideoFile:
    def __init__(self, line, directory):
        mouse_number, self.start_time, self.end_time, self.framerate, self.duration = line.split(",")
        self.base_path = os.path.join(directory, mouse_number)
        if os.path.exists(self.base_path + ".mkv"):
            self.video_file_path = self.base_path + ".mkv"
        elif os.path.exists(self.base_path + ".mp4"):
            self.video_file_path = self.base_path + ".mp4"
        else:
            raise IOError("No video file at %s" % self.base_path)
        self.still_path = self.base_path + " Still.jpg"
        self.mask_file_path = self.base_path + " Mask.jpg"
        self.start_time = datetime.datetime.strptime(self.start_time, "%m/%d/%Y %I:%M:%S %p")
        self.end_time = datetime.datetime.strptime(self.end_time, "%m/%d/%Y %I:%M:%S %p")
        self.movement_file_path = os.path.join("Movement Files", mouse_number + ".tsv")


class FrameToTime:
    """
    Designed to create a map between the frame number and the user input times
    """
    def __init__(self, start_time, end_time, framerate, duration):
        self.start_time = start_time
        self.end_time = end_time
        self.framerate = int(framerate)
        #t = datetime.datetime.strptime(duration, "%H:%M:%S")
        #self.duration = datetime.timedelta(hours=t.hour, minutes=t.minute, seconds=t.second)
        self.duration = pd.to_timedelta(duration)
        self.fps = self.framerate * self.duration.seconds / (self.end_time - self.start_time).total_seconds()
        print("FPS: %.2f" % self.fps)

    def set_number_of_frames(self, number_of_frames):
        self.fps = number_of_frames / (self.end_time - self.start_time).total_seconds()

    def frame_to_time(self, frame_number):
        time_after_start = datetime.timedelta(seconds=frame_number / self.fps)
        return self.start_time + time_after_start


def load_mask(mask_file_path):
    im = Image.open(mask_file_path)
    mask = np.array(im)
    mask = imutils.resize(mask, width=500)
    if len(mask.shape) > 2:
        mask = cv2.cvtColor(mask, cv2.COLOR_BGR2GRAY)
    mask = cv2.threshold(mask, 25, 255, cv2.THRESH_BINARY)[1]
    return mask


def create_sirenia_score_file(record, threshold, output_file):
    """
    This attempts to score each epoch as sleep or wake based off the number of movement pixels, and write that
    information to a file that looks like the score tsv file exported by Sirenia
    :param record: list of tuples, of the form (epoch start datetime.datetime, number of changed pixels in that epoch)
    :param threshold: number of changed pixels per epoch; epochs with fewer will be scored as sleep
    :param output_file: filepath to write to
    :return:
    """
    with open(output_file, "w+") as f:
        start_count = datetime.datetime(year=1899, month=12, day=30, hour=19, minute=0, second=0)
        f.write("Channels:\t1\n")
        f.write("Count:\t%d\n" % len(record))

        time_since_start_count = record[0][0] - start_count
        f.write("Start:\t%.8f\t%s\n" %
                (time_since_start_count.days + time_since_start_count.seconds / 86400,
                 record[0][0].strftime("%m/%d/%Y %I:%M:%S %p")))
        time_since_start_count = record[-1][0] - start_count
        f.write("End:\t%.8f\t%s\n" %
                (time_since_start_count.days + time_since_start_count.seconds / 86400,
                 record[-1][0].strftime("%m/%d/%Y %I:%M:%S %p")))
        f.write("Parameters\t4\n")
        f.write("NonRem\t2\n")
        f.write("REM\t3\n")
        f.write("Unscored\t255\n")
        f.write("Wake\t1\n\n")
        f.write("Date\tTime\tTime Stamp\tTime from Start\tCompy_Numeric\n")

        for i, r in enumerate(record):
            date = r[0].strftime("%m/%d/%Y")
            time = r[0].strftime("%H:%M:%S")
            time_since_start_count = r[0] - start_count
            time_stamp = time_since_start_count.days + time_since_start_count.seconds / 86400
            time_from_start = 10 * i
            score = "2" if r[1] < threshold else "1"
            f.write("%s\t%s\t%.8f\t%d\t%s\n" % (date, time, time_stamp, time_from_start, score))


def create_movement_file(record, output_file):
    """
    Instead of attempting to score each epoch as sleep or wake, this function just exports the information into a
    tsv file, listing the epoch and how many pixels changed during it
    :param record: list of tuples, of the form (epoch start datetime.datetime, number of changed pixels in that epoch)
    :param output_file: filepath to write to
    :return:
    """
    with open(output_file, "w+") as f:
        for r in record:
            f.write("%s\t%i\n" % (r[0].strftime("%m/%d/%Y %H:%M:%S"), r[1]))


def generate_still(video_file_path):
    vs = cv2.VideoCapture(video_file_path)
    last_frame = vs.read()[1]
    last_frame = imutils.resize(last_frame, width=500)
    last_frame = cv2.cvtColor(last_frame, cv2.COLOR_BGR2GRAY)
    im = Image.fromarray(last_frame)
    im.save("Still.jpg")


def save_video_still(vf):
    mouse_number = os.path.splitext(os.path.basename(vf.video_file_path))[0]
    vs = cv2.VideoCapture(vf.video_file_path)
    if vs is None or not vs.isOpened():
        raise IOError("%s does not exist" % vf.video_file_path)
    frame = vs.read()[1]
    frame = imutils.resize(frame, width=500)
    frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
    cv2.imwrite(vf.still_path, frame)


def analyze_video(vf):
    ftm = FrameToTime(vf.start_time, vf.end_time, vf.framerate, vf.duration)
    EPOCH_LENGTH = datetime.timedelta(seconds=10)
    mouse_number = os.path.splitext(os.path.basename(vf.video_file_path))[0]
    print(mouse_number)
    vs = cv2.VideoCapture(vf.video_file_path)
    last_frame = vs.read()[1]
    last_frame = imutils.resize(last_frame, width=500)
    last_frame = cv2.cvtColor(last_frame, cv2.COLOR_BGR2GRAY)
    frame_number = 0
    mask = load_mask(vf.mask_file_path)

    """
    This section of code throws away video frames until we get to one that is lined up with
    the start of an EEG epoch
    if args["epoch"] is not None:
        epoch_start = datetime.datetime.strptime(args["epoch"], "%m/%d/%Y %I:%M:%S %p")
        while (ftm.frame_to_time(frame_number) - epoch_start).seconds % 10 != 0:
            frame = vs.read()[1]
            frame_number += 1
    """

    last_frame = np.multiply(last_frame, mask)
    record = []
    reached_end = False
    start = time.time()
    predicted_number_frames = ftm.framerate * ftm.duration.seconds
    start_record_time = ftm.frame_to_time(frame_number)
    epochs_from_start = 0
    while not reached_end:
        movement_pixels = 0
        epochs_from_start += 1
        frames = []
        b = np.repeat(mask[:, :, np.newaxis], 3, axis=2)
        while ftm.frame_to_time(frame_number) < start_record_time + EPOCH_LENGTH * epochs_from_start:
            text = "No movement"

            frame = vs.read()[1]
            if frame_number % 1000 == 0:
                number_done = int(30 * frame_number / predicted_number_frames)
                print("[%s%s] %s/%s" % ("#" * number_done, "-" * (30-number_done), frame_number, predicted_number_frames))
            if frame is None:
                reached_end = True
                break
            frame = imutils.resize(frame, width=500)
            frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

            #This code displays each frame individually to allow for frame-by-frame analysis


            frame = np.multiply(frame, mask)

            frame_delta = cv2.absdiff(last_frame, frame)
            thresh = cv2.threshold(frame_delta, 25, 255, cv2.THRESH_BINARY)[1]
            moved_pixels = np.sum(thresh) / 255
            movement_pixels += moved_pixels

            if frame_number % 1000 == 0 and moved_pixels >= 1:
                fig, ax = plt.subplots(nrows=2, ncols=2)
                fig.suptitle(start_record_time + EPOCH_LENGTH * (epochs_from_start - 1))
                ax[0, 0].imshow(last_frame, cmap="gray")
                ax[0, 1].imshow(frame, cmap="gray")
                ax[1, 0].imshow(frame_delta, cmap="gray")
                ax[1, 1].imshow(thresh, cmap="gray")
                plt.show()

            frame_number += 1
            last_frame = frame

        """
        This code outputs each epoch, along with the number of movement pixels
        if epochs_from_start % 30 == 0:
            fourcc = cv2.VideoWriter_fourcc(*'DIVX')
            size = frames[0].shape[1], frames[0].shape[0]  ## <<<--- NOTICE THIS
            out = cv2.VideoWriter('Epochs/%s %06.f.avi' % (mouse_number, movement_pixels / 255)
                                  , fourcc, 12.0, size)
            for frame in frames:
                out.write(frame)
            out.release()
        """
        record.append((frame_number, movement_pixels))
    print("Frame Number", frame_number)
    print("Predicted Frame Number", predicted_number_frames)
    ftm.set_number_of_frames(frame_number)
    modified_record = []
    for frame_number, movement_pixels in record:
        modified_record.append((ftm.frame_to_time(frame_number), movement_pixels))
    create_movement_file(modified_record, vf.movement_file_path)
    print("Time: %.2f s" % (time.time() - start))


if __name__ == "__main__":
    """
    This program analyzes videos of mouse movement, outputting the amount of movement to a tsv file
    The metadata file should have lines in the following format:
    Mouse number, Start Time, End Time, Framerate, Video Duration
    Example:
    02550,9/27/2019  5:51:01 AM,9/29/2019  11:22:20 AM,12,13:24:04
    
    The video files should be saved as Mouse number.mkv or Mouse number.mp4 (Ex: 02550.mkv)
    
    This program will create save stills from each video file, saved as Mouse number Still.jpg (Ex: 02550 Still.jpg)
    You will need to open these images with photoshop, and create a black and white mask of the same dimensions.
    Black corresponds to areas that should not be analyzed, such as other cages, and white is areas that should
    be analyzed
    These masks should be saved as Mouse number Mask.jpg (Ex: 02550 Mask.jpg)
    The program will pause, waiting for input. Once you have finished creating the masks, press enter, and the program
    will continue creating the csv files.
    
    The csv files will be saved in Movement Files/Mouse number.tsv (Ex: Movement Files/02550.tsv). 
    Each line will have a format of
    Datetime\tMovement pixels
    (Ex: 02/22/2020 07:00:00 AM,536)
    """
    ap = argparse.ArgumentParser()
    ap.add_argument("-m", "--metadata", help="Path to metadata file")
    args = vars(ap.parse_args())
    videos = []
    directory = os.path.dirname(args["metadata"])
    with open(args["metadata"]) as f:
        for line in f:
            videos.append(VideoFile(line.strip(), directory))
    for vf in videos:
        save_video_still(vf)
    _ = input("After you have finished masking, press enter to continue: ")
    for vf in videos:
        if not os.path.exists(vf.movement_file_path):
            analyze_video(vf)