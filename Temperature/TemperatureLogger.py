import math
import serial
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import datetime
import requests
from slack_bolt import App
from slack_bolt.adapter.socket_mode import SocketModeHandler


def send_slack_message(message):
    dict = {'token': "xoxb-3822418182772-3832704777505-FIAU1YqRKdzGABHwnB57e9GN",
                   'channel': "#mouse-temperature-logging",
                   'text': message,
                   "icon_emoji": ":see_no_evil:",
                   "username": "Wobulator",
                   "blocks": None}
    requests.post('https://slack.com/api/chat.postMessage', dict).json()
    print(message)

class EMitter:
    def __init__(self, emitter_id):
        self.emitter_id = emitter_id
        self.m = 1
        self.b = 0
        if self.emitter_id == "EM00279":
            self.emitter_a = (.25)*np.log(984.3/822.5)
            self.emitter_b = np.log(984.3)-((41/4)*np.log(984.3/822.5))
        elif self.emitter_id == "EM01127":
            self.emitter_a = (1/(40.10-36.65))*np.log(947.1/811.78)
            self.emitter_b = np.log(947.1)-((40.10/(40.10-36.65))*np.log(947.1/811.78))
        elif self.emitter_id == "EM01328":
            self.emitter_a = (1/(41.0-36.6))*np.log(980.34/806.22)
            self.emitter_b = np.log(980.34)-((41.0/(41.0-36.6))*np.log(980.34/806.22))
        elif self.emitter_id == "EM01581":
            self.emitter_a = (1/(41-37))*np.log(977.6/817.2)
            self.emitter_b = np.log(977.6)-((41/(41-37))*np.log(977.6/817.2))
        elif self.emitter_id == "EM01584":
            self.emitter_a = (1/(41-37))*np.log(999.1/835.8)
            self.emitter_b = np.log(999.1)-((41/(41-37))*np.log(999.1/835.8))
        elif self.emitter_id == "EM01586":
            self.emitter_a = (1/(41-37))*np.log(962.4/804.5)
            self.emitter_b = np.log(962.4)-((41/(41-37))*np.log(962.4/804.5))
        elif self.emitter_id == "EM01833":
            self.emitter_a = (1/(41-37))*np.log(979.9/819.7)
            self.emitter_b = np.log(979.9)-((41/(41-37))*np.log(979.9/819.7))
        elif self.emitter_id == "EM01834":
            self.emitter_a = (1/(41-37))*np.log(964.8/806.7)
            self.emitter_b = np.log(964.8)-((41/(41-37))*np.log(964.8/806.7))
        elif self.emitter_id == "EM01836":
            self.emitter_a = (1/(41-37))*np.log(981.6/822.1)
            self.emitter_b = np.log(981.6)-((41/(41-37))*np.log(981.6/822.1))
        elif self.emitter_id == "EM01837":
            self.emitter_a = (1/(41-37))*np.log(1010.7/845)
            self.emitter_b = np.log(1010.7)-((41/(41-37))*np.log(1010.7/845))
        elif self.emitter_id == "EM03625":
            self.emitter_a = (1/(41-37))*np.log(985.1/823)
            self.emitter_b = np.log(985.1)-((41/(41-37))*np.log(985.1/823))
        elif self.emitter_id == "EM03626":
            self.emitter_a = (1/(41-37))*np.log(994.9/831.7)
            self.emitter_b = np.log(994.9)-((41/(41-37))*np.log(994.9/831.7))
        elif self.emitter_id == "EM03627":
            self.emitter_a = (1/(41-37))*np.log(866.6/724.6)
            self.emitter_b = np.log(866.6)-((41/(41-37))*np.log(866.6/724.6))
        elif self.emitter_id == "EM03628":
            self.emitter_a = (1/(41-37))*np.log(997.6/833.3)
            self.emitter_b = np.log(997.6)-((41/(41-37))*np.log(997.6/833.3))
        elif self.emitter_id == "EM03630":
            self.emitter_a = (1/(41-37))*np.log(963/805)
            self.emitter_b = np.log(963)-((41/(41-37))*np.log(963/805))
        elif self.emitter_id == "EM03656":
            self.emitter_a = (1/(41-37))*np.log(998.9/835.3)
            self.emitter_b = np.log(998.9)-((41/(41-37))*np.log(998.9/835.3))
        elif self.emitter_id == "EM03657":
            self.emitter_a = (1/(41-37))*np.log(986.1/823.5)
            self.emitter_b = np.log(986.1)-((41/(41-37))*np.log(986.1/823.5))
        elif self.emitter_id == "EM03658":
            self.emitter_a = (1/(41-37))*np.log(963.1/805.6)
            self.emitter_b = np.log(963.1)-((41/(41-37))*np.log(963.1/805.6))
        elif self.emitter_id == "EM03659":
            self.emitter_a = (1/(41-37))*np.log(982.8/821.4)
            self.emitter_b = np.log(982.8)-((41/(41-37))*np.log(982.8/821.4))
        elif self.emitter_id == "EM03660":
            self.emitter_a = (1/(41-37))*np.log(974/813.6)
            self.emitter_b = np.log(974)-((41/(41-37))*np.log(974/813.6))
            self.m = 1.0303
            self.b = -2.2523
        elif self.emitter_id == "EM03661":
            self.emitter_a = (1/(41-37))*np.log(984.8/823.4)
            self.emitter_b = np.log(984.8)-((41/(41-37))*np.log(984.8/823.4))
            self.m = 0.9857
            self.b = -0.185
        elif self.emitter_id == "EM03662":
            self.emitter_a = (1/(41-37))*np.log(1036.9/867.2)
            self.emitter_b = np.log(1036.9)-((41/(41-37))*np.log(1036.9/867.2))
            self.m = 1.0729
            self.b = -2.7435
        elif self.emitter_id == "EM03663":
            self.emitter_a = (1/(41-37))*np.log(976.8/816.2)
            self.emitter_b = np.log(976.8)-((41/(41-37))*np.log(976.8/816.2))
            self.m = 1.0074
            self.b = -0.5788
        elif self.emitter_id == "EM03664":
            self.emitter_a = (1/(41-37))*np.log(987.8/826.2)
            self.emitter_b = np.log(987.8)-((41/(41-37))*np.log(987.8/826.2))
        elif self.emitter_id == "EM03665":
            self.emitter_a = (1/(41-37))*np.log(967.9/809.6)
            self.emitter_b = np.log(967.9)-((41/(41-37))*np.log(967.9/809.6))
        elif self.emitter_id == "EM04182":
            self.emitter_a = (1/(41-37))*np.log(989.6/826.6)
            self.emitter_b = np.log(989.6)-((41/(41-37))*np.log(989.6/826.6))
            self.m = 1.0175
            self.b = -1.705
        elif self.emitter_id == "EM04183":
            self.emitter_a = (1/(41-37))*np.log(976.7/816.3)
            self.emitter_b = np.log(976.7)-((41/(41-37))*np.log(976.7/816.3))            
        elif self.emitter_id == "EM04184":
            self.emitter_a = (1/(41-37))*np.log(966.9/809.2)
            self.emitter_b = np.log(966.9)-((41/(41-37))*np.log(966.9/809.2))
        elif self.emitter_id == "EM04185":
            self.emitter_a = (1/(41-37))*np.log(1021.5/854.3)
            self.emitter_b = np.log(1021.5)-((41/(41-37))*np.log(1021.5/854.3))
            self.m = 1.0336
            self.b = -1.7123
        elif self.emitter_id == "EM04186":
            self.emitter_a = (1/(41-37))*np.log(984/822.9)
            self.emitter_b = np.log(984)-((41/(41-37))*np.log(984/822.9))
        elif self.emitter_id == "EM04187":
            self.emitter_a = (1/(41-37))*np.log(992.2/830)
            self.emitter_b = np.log(992.2)-((41/(41-37))*np.log(992.2/830))
            self.m = 0.907
            self.b = 2.5925
        elif self.emitter_id == "EM04188":
            self.emitter_a = (1/(41-37))*np.log(1034.6/865.8)
            self.emitter_b = np.log(1034.6)-((41/(41-37))*np.log(1034.6/865.8))
            self.m = 1.112
            self.b = -4.8304
        elif self.emitter_id == "EM04189":
            self.emitter_a = (1/(41-37))*np.log(982.4/821.7)
            self.emitter_b = np.log(982.4)-((41/(41-37))*np.log(982.4/821.7))
            self.m = 0.9781
            self.b = 0.0323
        elif self.emitter_id == "EM04190":
            self.emitter_a = (1/(41-37))*np.log(959.7/803.3)
            self.emitter_b = np.log(959.7)-((41/(41-37))*np.log(959.7/803.3))
            self.m = 0.9901
            self.b = 0.0481
        elif self.emitter_id == "EM04191":
            self.emitter_a = (1/(41-37))*np.log(963.9/806.8)
            self.emitter_b = np.log(963.9)-((41/(41-37))*np.log(963.9/806.8))
            self.m = 1.0435
            self.b = -1.7962
        elif self.emitter_id == "EM04919":
            self.emitter_a = (1/(41-37))*np.log(1009.7 / 844.2)
            self.emitter_b = np.log(1009.7)-((41/(41-37))*np.log(1009.7 / 844.2))
            self.m = 1.0522
            self.b = -2.4825
        elif self.emitter_id == "EM04920":
            self.emitter_a = (1 / (41 - 37)) * np.log(982.1 / 821.2)
            self.emitter_b = np.log(982.1) - ((41 / (41 - 37)) * np.log(982.1 / 821.2))
            self.m = 0.9774
            self.b = 0.2831
        else:
            raise IOError("%s is not a valid EMitter" % self.emitter_id)

    def temperature(self, bytes):
        return (1 / self.emitter_a) * (np.log(9825000 / ((256 * bytes[1]) + bytes[2])) - self.emitter_b)


class Mouse:
    """
    One mouse can move between different ER4000s, but with the same E-mitter, output file, and axis
    """
    def __init__(self, Emitter, ER4000s):
        if type(ER4000s) == list:
            self.ER4000s = ER4000s
        else:
            self.ER4000s = [ER4000s]
        self.Emitter = Emitter
        self.output_file = open("%s.csv" % self.__str__(), "a+")
        self.ax = None
        self.last_temperatures = []

    def __str__(self):
        return "-".join(self.ER4000s)

    def get_temperature(self):
        valid_temperatures = []
        for temperature in reversed(self.last_temperatures):
            if temperature != 0:
                valid_temperatures.append(temperature)
            if len(valid_temperatures) > 10:
                break
        if len(valid_temperatures) == 0:
            return 0
        else:
            return sum(valid_temperatures) / len(valid_temperatures)


def get_line(ser):
    return ser.readline().decode("utf-8").strip()


def close(event):
    send_slack_message("Logger stopped")

class Animator:
    def __init__(self, ser, mice):
        self.ser = ser
        self.mice = mice
        self.last_notification_time = 0
        """
        Different ER-4000 units return slightly different bytes. As far as I can figure out, ER-4000 units purchased
        from Ebay with serial codes that end in "00" send back 9 bytes, whereas the other units from the Buchanan
        lab send back 10 bytes. The extra byte is a checksum
        self.byte_length_lookup is a dictionary, that matches the ER-4000 to the number of bytes sent back
        """
        self.byte_length_lookup = {}
        for i in range(20):
            if i == 15 or i == 16:
                byte_length = 10
            else:
                byte_length = 9
            self.byte_length_lookup[str(i)] = byte_length
        self.total_ER4000s = 0
        for mouse in self.mice:
            self.total_ER4000s += len(mouse.ER4000s)
        self.ser.write(b"\t\r\n")
        bytes = self.ser.read(10 * self.total_ER4000s)
        self.last_interval = None

    def animate(self, _):
        self.last_interval = time.time()
        self.ser.write(b"\t\r\n")
        bytes = self.ser.read(10 * self.total_ER4000s)
        mice_temperatures_readings = {} # mice that already got a valid temperature reading this read
        for i in range(self.total_ER4000s):
            id = str(bytes[0])
            if id not in self.byte_length_lookup:
                return
            byte_block = bytes[:self.byte_length_lookup[id]]
            bytes = bytes[self.byte_length_lookup[id]:]
            mouse = None
            for iterate_mouse in self.mice:
                if id in iterate_mouse.ER4000s:
                    mouse = iterate_mouse
            else:
                pass
            if len(byte_block) == 10:  # If this is an E-mitter that sends a checksum
                if sum(byte_block[0:9]) % 256 != byte_block[9]:  # If the checksum is not right
                    return
            if byte_block[8] == 1 or byte_block[8] == 3:
                temperature = mouse.Emitter.temperature(byte_block)
                if mouse in mice_temperatures_readings.keys():
                    print("----------")
                    print("%s HAVE MULTIPLE READINGS" % str(mouse))
                    print("----------")
                else:
                    mice_temperatures_readings[mouse] = temperature
                    print(" ".join([str(b) for b in byte_block]))
            else:
                print("Invalid temperature")
        for mouse in self.mice:
            if mouse in mice_temperatures_readings:
                temperature = mice_temperatures_readings[mouse]
            else:
                temperature = 0
            mouse.last_temperatures.append(temperature)
            if len(mouse.last_temperatures) > 1000:
                mouse.last_temperatures.pop(0)
            if (0 < temperature < 20) or temperature > 42:
                if time.time() - self.last_notification_time > 5 * 60:
                    self.last_notification_time = time.time()
                    send_slack_message("Cage %s is at %.2f C" % (str(mouse), temperature))
            mouse.ax.clear()
            mouse.ax.set_title("%s: %.2f" % (str(mouse), temperature))
            mouse.ax.set_ylim(15, 45)
            mouse.ax.plot(np.arange(len(mouse.last_temperatures)), mouse.last_temperatures, 'r-')
            mouse.ax.set_ylabel("Temperature (C)")
            mouse.ax.set_xticks([])
            mouse.output_file.write("%s,%s\n" % (datetime.datetime.now(), temperature))


if __name__ == "__main__":
    bot_token = "xoxb-3822418182772-3832704777505-FIAU1YqRKdzGABHwnB57e9GN"
    app_token = "xapp-1-A03QGKNLT1P-3813461582710-e4f88b211058d038649921c77b228f556b6b5bcc601b7edb56f5a9a50d212ce1"
    app = App(token=bot_token)
    SocketModeHandler(app, app_token).connect()
    send_slack_message("Hello")

    mice = [
        Mouse(EMitter("EM04182"), ["1", "3", "5"]),
        Mouse(EMitter("EM04190"), ["2", "4", "6"])]
    ser = serial.Serial(port="COM5", baudrate=19200, timeout=1, bytesize=serial.EIGHTBITS)
    for mouse in mice:
        for ER4000 in mouse.ER4000s:
            ascii_character = chr(int(ER4000) + 48)
            ser.write(str.encode("%s,filt,4\r\n" % ascii_character))
            if get_line(ser) != "OK":
                raise IOError("%s did not respond to setting the filter" % ER4000)
            ser.write(str.encode("%s,hfilt,4\r\n" % ascii_character))
            if get_line(ser) != "OK":
                raise IOError("%s did not respond to setting the heart filter" % ER4000)
            ser.write(str.encode("%s,gain,2\r\n" % ascii_character))
            if get_line(ser) != "OK":
                raise IOError("%s dit not respond to setting the gain" % ER4000)
            ser.write(str.encode("%s,init\r\n" % ascii_character))
            time.sleep(2)
            if get_line(ser) != "OK":
                raise IOError("%s did not respond to initialization" % ER4000)
    send_slack_message("Logger started")
    last_temperatures = []
    plt.show()
    ncols = int(math.sqrt(len(mice)))
    nrows = int(math.sqrt(len(mice)) + 0.99)
    fig, axs = plt.subplots(ncols=ncols, nrows=nrows, squeeze=False)
    for mouse, ax in zip(mice, axs.reshape(-1)):
        mouse.ax = ax
    a = Animator(ser, mice)

    @app.message("Check")
    def message_hello(message, say):
        for mouse in a.mice:
            say("%s:%.2f" % (str(mouse), mouse.get_temperature()))

    _ = anim.FuncAnimation(fig=fig, func=a.animate, interval=1000, blit=False)
    plt.tight_layout()
    cid = fig.canvas.mpl_connect("close_event", close)
    plt.show()
    for mouse in mice:
        mouse.output_file.close()