import glob
import math
import os.path

import matplotlib.pyplot as plt
import matplotlib.patches
import numpy as np
import cv2
from PIL import Image
import pandas as pd
from matplotlib import font_manager
import skimage.draw
import scipy.ndimage.filters
from matplotlib.widgets import TextBox
import tkinter
from tkinter import filedialog

"""
This program is used to analyze the data from the "temperature bridge" thermal gradient experiment. This temperature
bridge is constructed by putting a sheet of copper between a Peltier cooler and a hot plate, so that the cold end
is at approximately 15 C, and the hot end is at 60 C. Mice are placed in lanes on this sheet, and are recorded with
a timelapse from a visual camera. A thermal image of the temperature bridge is taken, preferably at the beginning and end.

This program then finds the mouses' position using an intensity-based search, then
converts the mouse's position into the temperature on which they were sitting. Doing this requires the user's help
registering the two images together. Click on one image, then click on the corresponding point on the second image.
Registration will require at least four pairs of points for homographic (8D, location, size, rotation, skew).
Double click on the last point to initiate the transformation.
"""


class Click:
    def __init__(self, ax, flir, still):
        self.colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k', 'r', 'g', 'b', 'c', 'm', 'y', 'k']
        self.ax = ax
        self.flir = flir
        self.still = still
        self.ready_for_second_click = False
        self.last_click_on_top = None
        self.last_click_coordinates = None
        self.top_coordinate_pairs = []
        self.bottom_coordinate_pairs = []
        self.finished = False
        self.warped_flir = None

    def onclick(self, event):
        if event.inaxes == self.ax[0]:
            circle = plt.Circle((event.xdata, event.ydata), 4, color=self.colors[len(self.bottom_coordinate_pairs)])
            self.ax[0].add_patch(circle)
            if self.ready_for_second_click:
                if not self.last_click_on_top:
                    self.top_coordinate_pairs.append((event.xdata, event.ydata))
                    self.bottom_coordinate_pairs.append(self.last_click_coordinates)
                    self.ready_for_second_click = False
            else:
                self.ready_for_second_click = True
            self.last_click_coordinates = (event.xdata, event.ydata)
            self.last_click_on_top = True
        elif event.inaxes == self.ax[1]:
            circle = plt.Circle((event.xdata, event.ydata), 4, color=self.colors[len(self.bottom_coordinate_pairs)])
            self.ax[1].add_patch(circle)
            if self.ready_for_second_click:
                if self.last_click_on_top:
                    self.top_coordinate_pairs.append(self.last_click_coordinates)
                    self.bottom_coordinate_pairs.append((event.xdata, event.ydata))
                    self.ready_for_second_click = False
            else:
                self.ready_for_second_click = True
            self.last_click_coordinates = (event.xdata, event.ydata)
            self.last_click_on_top = False
        plt.draw()
        if event.dblclick:
            H, _ = cv2.findHomography(np.array(self.top_coordinate_pairs), np.array(self.bottom_coordinate_pairs))
            self.warped_flir = cv2.warpPerspective(self.flir, H, (self.still.shape[1], self.still.shape[0]), borderMode=cv2.BORDER_REPLICATE)
            ax[0].imshow(self.warped_flir)
            still_overlay = cv2.adaptiveThreshold(cv2.cvtColor(still, cv2.COLOR_BGR2GRAY), 5, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY, 11, 2)
            warped_flir_without_overlay = cv2.warpPerspective(self.flir, H, (self.still.shape[1], self.still.shape[0]))
            ax[1].imshow(still_overlay + warped_flir_without_overlay)
            self.top_coordinate_pairs = []
            self.bottom_coordinate_pairs = []
            self.finished = True
            plt.draw()


def get_still_from_video(video_file_path):
    vs = cv2.VideoCapture(video_file_path)
    first_frame = vs.read()[1]
    first_frame = cv2.cvtColor(first_frame, cv2.COLOR_RGB2BGR)
    return np.array(first_frame)


def get_flir_image(flir_file_path, barrel_distortion=-2.6e-4):
    if os.path.isdir(flir_file_path):
        return np.mean([get_flir_image(x) for x in glob.glob(os.path.join(flir_file_path, "*.csv"))], axis=0)
    flir = np.genfromtxt(flir_file_path, delimiter=",", skip_header=6)
    flir = cv2.copyMakeBorder(src=flir, top=15, bottom=15, left=15, right=15,
                                          borderType=cv2.BORDER_REPLICATE)
    width = flir.shape[1]
    height = flir.shape[0]

    distCoeff = np.zeros((4, 1), np.float64)

    k1 = barrel_distortion
    k2 = 0.0
    p1 = 0.0
    p2 = 0.0

    distCoeff[0, 0] = k1
    distCoeff[1, 0] = k2
    distCoeff[2, 0] = p1
    distCoeff[3, 0] = p2

    # assume unit matrix for camera
    cam = np.eye(3, dtype=np.float32)

    cam[0, 2] = width / 2.0  # define center x
    cam[1, 2] = height / 2.0  # define center y
    cam[0, 0] = 10.  # define focal length x
    cam[1, 1] = 10.  # define focal length y

    # here the undistortion will be computed
    return cv2.undistort(flir, cam, distCoeff)


def load_ethovision_tracks(ethovision_file_path):
    tracks = {}
    with pd.ExcelFile(ethovision_file_path) as xlsx:
        for sheet_name in xlsx.sheet_names:
            df = xlsx.parse(sheet_name, header=32, skiprows=[33])
            tracks[sheet_name] = np.array([df["X center"] + 960, np.array(df["Y center"]) + 540]).astype(int)
    return tracks


class Arena:
    def __init__(self, still, ax, textbox, axbox):
        self.shape = (still.shape[0], still.shape[1])
        self.masks = {}
        self.ax = ax
        self.vertices = []
        self.polygon_number = 0
        self.colors = ['r', 'g', 'b', 'k']
        self.textbox = textbox
        self.axbox = axbox
        self.axbox.set_visible(False)
        self.current_mask = None

    def onclick(self, event):
        if self.current_mask is None:
            if event.inaxes == self.ax:
                circle = plt.Circle((event.xdata, event.ydata), 4, color=self.colors[self.polygon_number])
                self.ax.add_artist(circle)
                if len(self.vertices) > 0:
                    last_line = plt.Line2D([self.vertices[-1][0], event.xdata], [self.vertices[-1][1], event.ydata], color=self.colors[self.polygon_number])
                    self.ax.add_line(last_line)
                self.vertices.append([event.xdata, event.ydata])
                if event.dblclick:
                    final_line = plt.Line2D([self.vertices[0][0], event.xdata], [self.vertices[0][1], event.ydata], color=self.colors[self.polygon_number])
                    self.ax.add_line(final_line)
                    self.polygon_number += 1
                    self.vertices = np.array(self.vertices)
                    rr, cc = skimage.draw.polygon(self.vertices[:, 1], self.vertices[:, 0], self.shape)
                    self.current_mask = np.zeros(self.shape, dtype=np.uint8)
                    self.current_mask[rr, cc] = 1
                    self.vertices = []
                    self.axbox.set_visible(True)
                plt.draw()

    def submit(self, text):
        if self.current_mask is not None:
            mouse_number = text
            self.masks[mouse_number] = self.current_mask
            self.current_mask = None
            self.axbox.set_visible(False)
            self.textbox.set_val("")


def draw_arenas(still):
    fig, ax = plt.subplots()
    axbox = plt.axes([0.3, 0.05, 0.6, 0.075])
    text_box = TextBox(axbox, 'Mouse number:', initial="")
    arena = Arena(still, ax, text_box, axbox)
    ax.imshow(still)
    ax.set_xticks([])
    ax.set_yticks([])
    cid = fig.canvas.mpl_connect("button_press_event", arena.onclick)
    text_box.on_submit(arena.submit)
    ax.set_title("Click to draw arenas, double click to close polygon, close to finish")
    plt.show()
    return arena.masks


def find_mouse_position(video_file_path, masks):
    if len(masks.items()) == 0:
        return None
    positions = {mouse_number: [] for mouse_number, mask in masks.items()}
    shape = list(masks.values())[0].shape
    array = np.zeros(shape=shape, dtype=np.uint8)
    vs = cv2.VideoCapture(video_file_path)
    for _ in range(100):
        frame = vs.read()[1]
        if frame is None:
            break
        frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        array = np.maximum(array, frame)
    plt.imshow(array, cmap="gray")
    plt.xticks([])
    plt.yticks([])
    plt.title("Acceptable empty arena?")
    plt.show()
    vs = cv2.VideoCapture(video_file_path)
    frame_number = 0
    while True:
        frame_number += 1
        frame = vs.read()[1]
        if frame is None:
            break
        frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        difference = cv2.absdiff(frame, array)
        thresh = cv2.threshold(difference, 75, 1, cv2.THRESH_BINARY)[1]
        minimum = cv2.bitwise_not(cv2.threshold(frame, 60, 255, cv2.THRESH_BINARY)[1])
        thresh = cv2.multiply(thresh, minimum)
        cXs, cYs = [], []
        for mouse_number, mask in masks.items():
            masked = cv2.multiply(thresh, mask)
            contours, hierarchy = cv2.findContours(masked, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
            if len(contours) > 0:
                c = max(contours, key=cv2.contourArea)
                M = cv2.moments(c)
                if M["m00"] == 0:
                    fig, ax = plt.subplots(3)
                    ax[0].imshow(frame, cmap="gray")
                    ax[1].imshow(difference, cmap="gray")
                    ax[2].imshow(thresh, cmap="gray")
                    plt.show()
                cX = int(M["m10"] / M["m00"])
                cY = int(M["m01"] / M["m00"])
                cXs.append(cX)
                cYs.append(cY)
                positions[mouse_number].append([cX, cY])
        if frame_number % 500 == 0:
            fig, ax = plt.subplots()
            ax.imshow(frame, cmap="gray")
            for cX, cY in zip(cXs, cYs):
                circ = matplotlib.patches.Circle((cX, cY), 4, facecolor='r')
                ax.add_patch(circ)
            ax.set_xticks([])
            ax.set_yticks([])
            plt.title(frame_number)
            plt.show()
    return {mouse_number: np.array(position) for mouse_number, position in positions.items()}


if __name__ == "__main__":
    font_path = r"C:\Users\Fillan\Documents\Fonts\Myriad Pro Regular.ttf"
    font_manager.fontManager.addfont(font_path)
    prop = font_manager.FontProperties(fname=font_path)
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = prop.get_name()
    plt.rcParams['font.size'] = 14
    plt.rcParams['axes.linewidth'] = 2

    root = tkinter.Tk()
    root.withdraw()
    folder = tkinter.filedialog.askdirectory()
    video_file_path = glob.glob(os.path.join(folder, "*.AVI"))[0]
    files = glob.glob(os.path.join(folder, "*.csv"))
    flir = get_flir_image(folder)
    display_flir = flir - scipy.ndimage.filters.gaussian_filter(flir, sigma=20)
    nrows = int(math.sqrt(len(files)))
    fig, axs = plt.subplots(nrows=nrows, ncols=math.ceil(len(files) / nrows))
    for ax, file in zip(axs.reshape(-1), files):
        ax.imshow(flir - get_flir_image(file))
        ax.set_xticks([])
        ax.set_yticks([])
    plt.show()
    still = get_still_from_video(video_file_path)
    arenas = draw_arenas(still)
    positions = find_mouse_position(video_file_path, arenas)
    fig, ax = plt.subplots(nrows=2)
    im = Image.fromarray(flir)
    im = im.convert("L")
    im.save(os.path.join(folder, "Flir.jpg"))
    ax[0].imshow(display_flir)
    ax[1].imshow(still)
    ax[0].set_title("Click fiducial on top image, then corresponding fiducial on bottom. Double click to finish")
    click = Click(ax, flir, still)
    for axis in ax:
        axis.set_xticks([])
        axis.set_yticks([])
    cid = fig.canvas.mpl_connect("button_press_event", click.onclick)
    plt.tight_layout()
    plt.show()
    if click.finished:
        mouse_temperatures = {}
        for mouse_number, track in positions.items():
            print(mouse_number)
            temperatures = []
            with open(os.path.join(os.path.dirname(video_file_path), "%s.tsv" % mouse_number), "w+") as f:
                for i in range(track.shape[0]):
                    temperature = click.warped_flir[track[i][1], track[i][0]]
                    temperatures.append(temperature)
                    f.write("%.2f\t%.2f\n" % (i / 6, temperature))
            plt.plot([i / 6 for i in range(len(temperatures))], temperatures, label=mouse_number)
            mouse_temperatures[mouse_number] = temperatures
        plt.ylabel("Temperature")
        plt.xlabel("Minutes")
        plt.legend()
        plt.show()
        for mouse_number, temperatures in mouse_temperatures.items():
            plt.hist(temperatures, label=mouse_number, histtype="step", bins=np.arange(10, 45))
            print("%s: Mean: %.2f Stdev: %.2f" % (mouse_number, np.mean(temperatures), np.std(temperatures)))
        plt.legend()
        plt.ylabel("N")
        plt.xlabel("Temperature (C)")
        plt.show()



