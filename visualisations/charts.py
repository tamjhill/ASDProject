import matplotlib.pyplot as plt

def plot_size(sizetotal, articletotal):
    avesize = sizetotal/articletotal
    x = range(articletotal)
    y1 = [xi*avesize for xi in x]
    y2 = [xi*avesize*0.0403 for xi in x]
    y3 = [xi*avesize*0.048 for xi in x]
    y4 = [xi*avesize*0.0583 for xi in x]

    plt.plot(x, y1, label='Projected size before compression')
    plt.plot(x, y2, label='Predicted size after direct compression')
    plt.plot(x, y2, label='Predicted size after adjacency compression')
    plt.plot(x, y2, label='Predicted size after tree+triple compression')
    plt.legend()
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.show()
    plt.savefig('compression.png')

#plot triples
plot_size(38295740, 140)
#plot mb
plot_size(478965, 140)