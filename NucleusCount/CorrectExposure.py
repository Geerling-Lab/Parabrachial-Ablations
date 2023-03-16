import PIL.Image
import numpy as np
import glob
import matplotlib.pyplot as plt
TARGET_EXPOSURE = 180
TARGET_CONTRAST = 30

if __name__ == "__main__":
    files = glob.glob(r"R:\NucleusCount\05131B\*\Image.png")
    for file in files:
        img = PIL.Image.open(file)
        width, height = img.size
        arr = np.array(img)
        fig, ax = plt.subplots(ncols=2, sharex=True, sharey=True)
        ax[0].imshow(arr)
        exposure = np.mean(arr)
        print(file)
        print("Previous Exposure: %.2f" % exposure)
        contrast = np.std(arr)
        print("Previous Contrast: %.2f" % contrast)
        arr = (arr - exposure) * (TARGET_CONTRAST / contrast) + TARGET_EXPOSURE
        arr[arr > 255] = 255
        arr[arr < 0] = 0
        arr = arr.astype(np.int16)
        exposure = np.mean(arr)
        print(file)
        print("After Exposure: %.2f" % exposure)
        contrast = np.std(arr)
        print("After Contrast: %.2f" % contrast)
        ax[1].imshow(arr)
        for a in ax:
            a.set_xticks([])
            a.set_yticks([])
        plt.tight_layout()
        plt.show()