import matplotlib.pyplot as plt

f = open("out.txt", "r")


x = []
y = []

n = int(f.readline())

for n in range(n):
    s = f.readline()
    [xx, yy] = list(map(float, s.split()))
    x.append(xx)
    y.append(yy)

plt.plot(x, y, 'ro')
plt.axis([-10, 10, -10, 10])
plt.show()