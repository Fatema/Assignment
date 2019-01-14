import argparse
from math import floor, sqrt
from random import random

parser = argparse.ArgumentParser()
parser.add_argument("nbody", help="number of random bodies to generate")
args = parser.parse_args()

params = "0 1e4"
xoffset = 0
yoffset = 0
x = 0
y = 0
z = 0
m = 39.948
n = int(args.nbody)
limit = floor(sqrt(n))

for i in range(1, n + 1):
    x = "{:0.1e}".format(xoffset) if xoffset != 0 else 0
    y = "{:0.1e}".format(yoffset) if yoffset != 0 else 0
    xsign = -1 if random() < 0.5 else 1
    ysign = -1 if random() < 0.5 else 1
    vx = 1e-12 * xsign if random() < 0.5 else 0
    vy = 1e-12 * ysign if random() < 0.5 else 0
    vz = 0
    xoffset = xoffset + 5e-9 if i % limit > 0 else 0
    yoffset = yoffset + 5e-9 if i % limit == 0 else yoffset
    params += " {x} {y} {z} {vx} {vy} {vz} {m}".format(x=x, y=y, z=z, m=m, vx=vx, vy=vy, vz=vz)

print params


