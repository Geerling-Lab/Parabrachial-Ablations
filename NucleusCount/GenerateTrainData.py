import argparse
import numpy as np
import PIL.Image
import random
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import os
from Constants import SIZE, RADIUS
import scipy.ndimage

"""
This file creates the training data necessary for training the network. Give it an image file and csv files
with the human labeled structures.
"""


class CellMap:
    def __init__(self, arr, greyscale_arr):
        self.arr = arr
        self.greyscale_arr = greyscale_arr
        self.cell_map = np.zeros(greyscale_arr.shape, dtype=bool)
        self.counter = 0


def data_augmentation(X, Y, rotate=True, adjust_brightness=False, repeat=1):
    """
    Allows for data augmentation to increase the training set
    """
    Xs = [X]
    if rotate:
        Xs.append(np.rot90(X))
        Xs.append(np.rot90(X, k=2))
        Xs.append(np.rot90(X, k=3))
    if adjust_brightness:
        for i in range(5):
            a = (X * (1 + (i - 3) * 0.1)).astype(X.dtype)
            Xs.extend([a, a[::-1, :], a[:, ::-1], a[::-1, ::-1]])
    Xs = Xs * repeat
    Ys = [Y] * len(Xs)
    return Xs, Ys


def generate_negative_data(arr, greyscale_arr, cell_map, number):
    """
    Creates negative training data for training the network
    Creates NUMBER (x,y) pairs
    x is an array of size SIZE*SIZE drawn at random from arr, weighted by the intensity of that pixel in arr
    If x overlaps too much with the user selected positive cells from cell_map, throw it out
    :param arr:
    :param cell_map:
    :param number:
    :return:
    """
    X = []
    Y = []
    probabilities = (255 - np.ravel(scipy.ndimage.gaussian_filter(greyscale_arr, 10)).astype(np.uint64))
    probabilities = probabilities ** 3
    probabilities = probabilities / probabilities.sum()
    last_progress_bar_length = 0
    while len(X) < number:
        progress_bar_length = int(80 * len(X) / number)
        if progress_bar_length > last_progress_bar_length:
            print("[" + "-" * progress_bar_length + " " * (80 - progress_bar_length) + "]")
            last_progress_bar_length = progress_bar_length
        start_locations = np.random.choice(greyscale_arr.size, size=100, p=probabilities)
        # Generates a weighted random vector 100 long, where the probability of any value is the lightness of that pixel
        xs, ys = np.unravel_index(start_locations, greyscale_arr.shape)
        for x, y in zip(xs, ys):
            if np.average(cell_map[int(x - RADIUS): int(x + RADIUS), int(y - RADIUS): int(y + RADIUS)]) < 0.90:
                sub_arr = arr[int(x - RADIUS): int(x + RADIUS), int(y - RADIUS): int(y + RADIUS), :]
                if sub_arr.size == SIZE ** 2 * 3:
                    Xs, Ys = data_augmentation(sub_arr, [0, 1], rotate=True, adjust_brightness=False)
                    X.extend(Xs)
                    Y.extend(Ys)
    return X, Y


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train NN")
    parser.add_argument("-c", "--csv", type=str, help="Filepath of input csv file")
    parser.add_argument("-t", "--train", type=str, help="Directory to output training data")
    args = parser.parse_args()
    cell_maps = {}
    Xs, Ys = [], []
    with open(args.csv) as f:
        for line in f:
            image_file, _, _, _, x, y, w, h = line.strip().split(",")
            x, y, w, h = int(x), int(y), int(w), int(h)
            if image_file not in cell_maps.keys():
                arr = np.array(PIL.Image.open(image_file)).astype(np.uint8)
                greyscale_arr = np.average(arr, axis=2)
                cell_maps[image_file] = CellMap(arr, greyscale_arr)
            subset = cell_maps[image_file].arr[x: x + w, y: y + h, :]
            if subset.shape != (SIZE, SIZE, 3):
                continue
            sub_X, sub_Y = data_augmentation(subset, [1, 0], adjust_brightness=False)
            cell_maps[image_file].cell_map[x: x + w, y: y + h] = True
            cell_maps[image_file].counter += len(sub_X)
            Xs.extend(sub_X)
            Ys.extend(sub_Y)
    for image_file, cell_map in cell_maps.items():
        number_unfilled_squares = np.size(cell_map.cell_map) - np.count_nonzero(cell_map.cell_map)
        sub_X, sub_Y = generate_negative_data(cell_map.arr, cell_map.greyscale_arr, cell_map.cell_map, number_unfilled_squares // 300)
        Xs.extend(sub_X)
        Ys.extend(sub_Y)
    X = np.array(Xs)
    Y = np.array(Ys)
    print("%s negative" % np.sum(Y[:, 1]))
    print("%s positive" % np.sum(Y[:, 0]))
    if not os.path.isdir(args.train):
        os.mkdir(args.train)
    save_file_path = os.path.join(args.train, "%s.npz" % os.path.basename(args.csv[:-4]))
    np.savez_compressed(save_file_path, X=X, Y=Y)