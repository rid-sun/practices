import matplotlib.pyplot as plt
import numpy as np

def plot_util(x_axis_data, y_axis_data, xlabel, ylabel, fig_name):
    for x, y in zip(x_axis_data, y_axis_data):
        plt.text(x, y+0.3, '%.2f' % y, horizontalalignment="center", verticalalignment="center", fontsize=7.5)# y_axis_data1加标签数据

    plt.plot(x_axis_data, y_axis_data, 'b*--', alpha=0.5, linewidth=1, label='vol')# 'bo-'表示蓝色实线，数据点实心原点标注
    ## plot中参数的含义分别是横轴值，纵轴值，线的形状（'s'方块,'o'实心圆点，'*'五角星   ...，颜色，透明度,线的宽度和标签 ，

    plt.legend()  # 显示上面的label
    plt.xlabel(xlabel) # x_label
    plt.ylabel(ylabel)# y_label
    
    # plt.ylim(-1,1)#仅设置y轴坐标范围
    
    # 自动保存图片到本地
    plt.savefig(fig_name + '.png', dpi = 500)
    
    # 要注意show放在savefig之后，因为默认show之后会重新打开白板，如果这时再进行save那么会是白板
    # plt.show()

if __name__ == '__main__':
    x_axis_data = [1,2,3,4,5,6,7] # x
    y_axis_data = [68,69,79,71,80,70,66] # y
    xlabel = 'time'
    ylabel = 'number'
    fig_name = 'fig0'
    plot_util(x_axis_data, y_axis_data, xlabel, ylabel, fig_name)