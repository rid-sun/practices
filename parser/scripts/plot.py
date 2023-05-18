import matplotlib.pyplot as plt
import numpy as np
from pylab import mpl
# 设置显示中文字体
mpl.rcParams["font.sans-serif"] = ["SimHei"]
# 设置正常显示符号
mpl.rcParams["axes.unicode_minus"] = False

def plot_util(x_axis_data, y_axis_data, xlabel, ylabel, title, fig_name):
    all_color = ['#DCDCDC', '#FFDAB9', '#FFE4E1', '#2F4F4F', '#6495ED', '#00BFFF', '#00FFFF', '#66CDAA', '#98FB98', '#DCDCDC', '#FFDAB9', '#FFE4E1', '#2F4F4F', '#6495ED', '#00BFFF', '#00FFFF', '#66CDAA', '#98FB98']
    max_m = float('-inf')
    min_l = float('inf')
    fig, ax = plt.subplots()
    # print(len(y_axis_data))

    for i in range(len(y_axis_data)):
        max_m = max(max(y_axis_data[i]), max_m)
        min_l = min(min(y_axis_data[i]), min_l)
        ax.plot(x_axis_data, y_axis_data[i], all_color[i], label = 'x' + str(i + 1))# y_axis_data1加标签数据
    
    ax.legend(loc = 'upper left')  # 显示上面的label
    ax.set_xlabel(xlabel) # x_label
    ax.set_ylabel(ylabel) # y_label
    ax.set_title(title)  # 设置标题
    ax.set_ylim(min_l - 0.5, max_m + 0.5)
    
    # 自动保存图片到本地
    plt.savefig(fig_name + '.png', dpi = 300)

    # 要注意show放在savefig之后，因为默认show之后会重新打开白板，如果这时再进行save那么会是白板
    # plt.show()

if __name__ == '__main__':
    x_axis_data = [1,2,3,4,5,6,7] # x
    y_axis_data = [[68,69,79,71,80,70,66]] # y
    xlabel = 'time'
    ylabel = 'number'
    fig_name = 'fig0'
    title = 'test'
    plot_util(x_axis_data, y_axis_data, xlabel, ylabel, title, fig_name)