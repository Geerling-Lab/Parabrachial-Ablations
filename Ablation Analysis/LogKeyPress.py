import pyautogui
import keyboard
import pytesseract
import os
from PIL import Image, ImageOps, ImageFilter
import numpy as np

class KP:
    def __init__(self, f):
        self.f = f

    def key_press_small(self):
        image = pyautogui.screenshot(region=(1145, 2020, 100, 22))
        thresh = 200
        fn = lambda x: 255 if x > thresh else x
        im2 = image.resize((300, 72), resample=Image.NEAREST)
        im2.save(r"C:\Users\fgrady\Downloads\im.png")
        arr = np.array(im2)
        arr[arr.sum(axis=2) > 600, :] = 255
        arr2 = np.zeros(shape=(200, 600, 3), dtype=np.uint8)
        arr2[:, :, :] = 255
        arr2[64:136, 150:450, :] = arr
        im3 = Image.fromarray(arr2)
        text = pytesseract.image_to_string(image, config="--dpi 300 --psm 8")
        text = ''.join(c for c in text if c.isdigit()).strip()
        text2 = pytesseract.image_to_string(im2)
        text2 = ''.join(c for c in text2 if c.isdigit()).strip()
        print("Text: %s    Text2: %s" % (text, text2))
        self.f.write("%s,%s%s" % (text, "0", os.linesep))

    def key_press_large(self):
        image = pyautogui.screenshot(region=(450, 1758, 3300, 50))
        im2 = ImageOps.invert(image)
        im2.show()
        text = pytesseract.image_to_string(im2)
        print("Text: %s" % text)
        text = ''.join(c for c in text if c.isdigit())
        self.f.write("%s,%s%s" % (text, "1", os.linesep))


if __name__ == "__main__":
    output_csv = r"C:\Users\fgrady\Downloads\MicturitionScores.csv"
    f = open(output_csv, "w+")
    myKP = KP(f)
    keyboard.add_hotkey('k', myKP.key_press_small, args=())
    keyboard.add_hotkey('l', myKP.key_press_large, args=())
    keyboard.wait('esc')
    f.close()