import os
import glob


def get_size(path):
    """
    path should be a csv file, created from ResearchIR by exporting an ROI
    This function simply returns the number of nonzero entries in the file
    """
    count = 0
    with open(path) as f:
        for line in f:
            words = line.strip().split(",")
            for word in words:
                try:
                    if len(word) == 0:
                        continue
                    elif float(word) > 0:
                        count += 1
                except ValueError as ve:
                    print(word)
    return count


if __name__ == "__main__":
    CAGE_SIZE = 18 * 35  # Known size of cage bottoms, in cm^2
    folder_location = r"C:\Users\fgrady\Downloads"
    output_file = r"C:\Users\fgrady\Downloads\Output.csv"
    if os.path.exists(output_file):
        os.remove(output_file)
    os.chdir(folder_location)
    result = glob.glob("*.{}".format("csv"))
    box_size_pixels = None
    output = {}
    for file in result:
        path = os.path.join(folder_location, file)
        if "box" in file.lower():
            box_size_pixels = get_size(path)
    for file in result:
        path = os.path.join(folder_location, file)
        if "box" not in file.lower():
            number = os.path.splitext(file)[0]
            size = CAGE_SIZE * get_size(path) / box_size_pixels
            output[number] = size
    print("Box Size: %.2f" % box_size_pixels)
    print("Pixel Size: %.2f" % (CAGE_SIZE / box_size_pixels))
    with open(output_file, "w") as f:
        for key, value in output.items():
            f.write("%s,%.2f\n" % (key, value))
            print("%s,%.2f\n" % (key, value))