import datetime

"""
If a file has been scored in 10 second epochs but we need 5 second epochs
"""

if __name__ == "__main__":
    input_file_path = r"R:\EEG Scores\Optogenetic EEG Scores\Scores\05960 scores.tsv"
    output_file_path = r"R:\EEG Scores\Optogenetic EEG Scores\Scores\05960 scores 5 seconds.tsv"
    with open(input_file_path) as f:
        with open(output_file_path, "w+") as of:
            for line in f:
                words = line.split("\t")
                if len(words) == 6:
                    try:
                        date_time_str = words[0] + " " + words[1]
                        date_time = datetime.datetime.strptime(date_time_str, "%m/%d/%Y %H:%M:%S")
                        time_stamp = float(words[2].strip())
                        time_from_start = float(words[3].strip())
                        of.write("\t".join(words))
                        date_time += datetime.timedelta(seconds=5)
                        words[0] = datetime.datetime.strftime(date_time, "%m/%d/%Y")
                        words[1] = "     " + datetime.datetime.strftime(date_time, "%H:%M:%S")
                        words[2] = "       %.9f" % (time_stamp + 0.000057871)
                        words[3] = "%.6f" % (time_from_start + 5)
                        of.write("\t".join(words))
                    except ValueError as ve:
                        of.write(line)
                else:
                    of.write(line)
