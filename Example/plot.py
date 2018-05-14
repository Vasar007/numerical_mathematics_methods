import argparse

import matplotlib.pyplot as plt


def main():
    # Construct the argument parser and parse the arguments.
    ap = argparse.ArgumentParser()
    ap.add_argument("-f", "--file", required=True, help="data file name")
    args = vars(ap.parse_args())

    title = ""
    x_label = ""
    y_label = ""
    for data in args["file"].split(","):
        data_x1 = []
        data_x2 = []
        data_y1 = []
        data_y2 = []
        data_y3 = []
        data_y4 = []
        data_y5 = []
        data_y6 = []
        data_y7 = []
        data_y8 = []
        mode = ""
        params = ""
        with open(data, "r") as f:
            (mode, params, title, x_label, y_label) = next(f).split("|")
            if mode == "1p":
                for line in f:
                    (x, y) = line.split()
                    data_x1.append(float(x))
                    data_y1.append(float(y))
                if params == "def":
                    plt.plot(data_x1, data_y1)
                else:
                    plt.plot(data_x1, data_y1, params)
            elif mode == "2p":
                flag = False
                for line in f:
                    if line == "#\n":
                        flag = True
                        continue
                    (x, y) = line.split()
                    if flag:
                        data_x1.append(float(x))
                        data_y1.append(float(y))
                    else:
                        data_x2.append(float(x))
                        data_y2.append(float(y))
                if params == "def":
                    plt.plot(data_x1, data_y1, data_x2, data_y2)
                else:
                    parsed_params = params.split(";")
                    plt.plot(data_x1, data_y1, parsed_params[0],
                             data_x2, data_y2, parsed_params[1])
            elif mode == "4p":
                for line in f:
                    (x, y1, y2, y3, y4) = line.split()
                    data_x1.append(float(x))
                    data_y1.append(float(y1))
                    data_y2.append(float(y2))
                    data_y3.append(float(y3))
                    data_y4.append(float(y4))
                if params == "def":
                    plt.plot(data_x1, data_y1, data_x1, data_y2, data_x1,
                             data_y3, data_x1, data_y4)
                else:
                    parsed_params = params.split(";")
                    plt.plot(data_x1, data_y1, parsed_params[0],
                             data_x1, data_y2, parsed_params[1],
                             data_x1, data_y3, parsed_params[2],
                             data_x1, data_y4, parsed_params[3])
            elif mode == "8p":
                flag = False
                for line in f:
                    if line == "#\n":
                        flag = True
                        continue
                    (x, y1, y2, y3, y4) = line.split()
                    if flag:
                        data_x1.append(float(x))
                        data_y1.append(float(y1))
                        data_y2.append(float(y2))
                        data_y3.append(float(y3))
                        data_y4.append(float(y4))
                    else:
                        data_x2.append(float(x))
                        data_y5.append(float(y1))
                        data_y6.append(float(y2))
                        data_y7.append(float(y3))
                        data_y8.append(float(y4))
                if params == "def":
                    plt.plot(data_x1, data_y1, data_x1, data_y2, data_x1,
                             data_y3, data_x1, data_y4, data_x1, data_y5,
                             data_x1, data_y6, data_x1, data_y7, data_x1,
                             data_y8)
                else:
                    parsed_params = params.split(";")
                    plt.plot(data_x1, data_y1, parsed_params[0],
                             data_x1, data_y2, parsed_params[1],
                             data_x1, data_y3, parsed_params[2],
                             data_x1, data_y4, parsed_params[3],
                             data_x2, data_y5, parsed_params[4],
                             data_x2, data_y6, parsed_params[5],
                             data_x2, data_y7, parsed_params[6],
                             data_x2, data_y8, parsed_params[7],)

        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        plt.grid(True)
        plt.axis("auto")
        plt.show()


if __name__ == '__main__':
    main()
