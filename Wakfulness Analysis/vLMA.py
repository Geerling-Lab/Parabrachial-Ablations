import cv2
import numpy as np
import matplotlib.pyplot as plt
import datetime
import time
from matplotlib.widgets import TextBox
import threading


class Click:
    def __init__(self, axs, axbox, text_box):
        self.axs = axs
        self.click_values = {}
        self.axbox = axbox
        self.just_selected = None
        self.text_box = text_box

    def click(self, event):
        for i, ax in enumerate(self.axs):
            if event.inaxes == ax:
                self.just_selected = i
                self.axbox.set_visible(True)
                plt.draw()
                if event.dblclick:
                    assert len(self.click_values) > 0
                    plt.close()

    def submit(self, text):
        if self.just_selected is not None:
            self.axbox.set_visible(False)
            self.click_values[self.just_selected] = text
            self.axs[self.just_selected].set_visible(False)
            self.just_selected = None
            self.text_box.set_val("")
            plt.draw()


class Drawer():
    def __init__(self, ax, camera_ax, cap, file, camera_name):
        self.ax = ax
        self.camera_ax = camera_ax
        self.cap = cap
        self.movement = []
        self.last_frame = None
        self.file = file
        self.epoch_movements = []
        self.last_time = None
        self.camera_name = camera_name
        self.stopped = False
        self.start_time = time.time()

    def run(self):
        if not self.stopped:
            threading.Timer(0.333, self.run).start()
        if self.last_time is not None:
            print("Lag: %.2f" % (time.time() - self.last_time))
        self.last_time = time.time()
        ret, frame = self.cap.read()
        if not ret:
            print("Invalid frame")
            return
        frame = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        if self.last_frame is not None:
            frame_delta = cv2.absdiff(self.last_frame, frame)
            thresh = cv2.threshold(frame_delta, 25, 255, cv2.THRESH_BINARY)[1]
            moved_pixels = np.sum(thresh) / 255
            print("%s got frame" % self.camera_name)
            self.movement.append(moved_pixels)
            if len(self.movement) > 30:
                epoch_movement = sum(self.movement)
                self.movement = []
                self.file.write("%s,%i\n" % (datetime.datetime.now().strftime("%m/%d/%y %H:%M:%S"), epoch_movement))
                self.file.flush()
                self.epoch_movements.append(epoch_movement)
                time_since_start = (time.time() - self.start_time) / 60
                xs = np.arange(0, time_since_start, time_since_start / len(self.epoch_movements))
                self.ax.plot(xs, self.epoch_movements, "r-")
                self.ax.set_ylabel(self.camera_name)
                self.camera_ax.imshow(frame, cmap="Greys_r")
                if len(self.epoch_movements) % 10 == 0:
                    self.ax.clear()
                    self.camera_ax.clear()
        self.last_frame = frame
        plt.draw()

    def close(self, event):
        print("Closing")
        self.file.close()
        self.stopped = True


if __name__ == "__main__":
    cameras = {}
    for camera_idx in range(6):
        cap = cv2.VideoCapture(camera_idx)
        print(camera_idx)
        try:
            if cap.isOpened():
                ret, frame = cap.read()
                cap.release()
                frame = cv2.cvtColor(frame, cv2.COLOR_RGB2BGR)
                cameras[camera_idx] = frame
        except cv2.error:
            print("%i camera invalid" % camera_idx)
            cameras[camera_idx] = np.zeros(shape=(100, 100, 3))
            cap.release()
    fig, axs = plt.subplots(ncols=len(cameras))
    fig.canvas.mpl_disconnect(fig.canvas.manager.key_press_handler_id)
    for (camera_index, frame), ax in zip(cameras.items(), axs):
        ax.imshow(frame)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(camera_index)
    axbox = plt.axes([0.3, 0.05, 0.6, 0.075])
    text_box = TextBox(axbox, 'Camera name:', initial="")
    axbox.set_visible(False)
    click = Click(axs, axbox, text_box)
    cid = fig.canvas.mpl_connect("button_press_event", click.click)
    plt.suptitle("Click on the correct camera")
    text_box.on_submit(click.submit)
    plt.show()
    fig, axs = plt.subplots(ncols=len(click.click_values), nrows=2, squeeze=False, gridspec_kw={'width_ratios': [3, 2]})
    plt.subplots_adjust(wspace=0, hspace=0)
    axs[1, 0].set_xlabel("Minutes")
    for i, (camera_number, camera_name) in enumerate(click.click_values.items()):
        print("Camera *%s* is called *%s**" % (camera_number, camera_name))
        file = open("%s %s.txt" % (camera_name, datetime.datetime.now().strftime("%b %d %Y")), "a+")
        cap = cv2.VideoCapture(camera_number)
        while not cap.isOpened():
            print("%s not opened" % camera_number)
            cap.release()
            cap = cv2.VideoCapture(camera_number, cv2.CAP_DSHOW)
        cap.set(cv2.CAP_PROP_FRAME_WIDTH, 320)
        cap.set(cv2.CAP_PROP_FRAME_HEIGHT, 180)
        drawer = Drawer(axs[i, 0], axs[i, 1], cap, file, camera_name)
        axs[i, 1].set_xticks([])
        axs[i, 1].set_yticks([])
        axs[i, 1].set_xticklabels([])
        axs[i, 1].set_yticklabels([])
        cid = fig.canvas.mpl_connect("close_event", drawer.close)
        drawer.run()
    plt.tight_layout()
    plt.show()