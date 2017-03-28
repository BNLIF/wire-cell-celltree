#!/usr/bin/env python
'''
Make some plots using wire-cell celltree geometry files
'''
# channel: channel ID
# plane: plane number (0 == U, 1 == V, 2 == Z)
# wire: wire ID on that plane
# sx, sy, sz: starting (x,y,z) position of the wire in cm
# ex, ey, ez: ending (x,y,z) position of the wire in cm

from collections import namedtuple
import numpy
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


Wire = namedtuple("Wire","ch plane wire beg end")

class Wires(object):
    def __init__(self, filename = None):
        if filename:
            self.load(filename)

    def load(self, filename):
        wires = list()
        with open(filename) as fp:
            for line in fp.readlines():
                if line.startswith("#"):
                    continue
                line = line.strip()
                if not line:
                    continue
                chunks = line.split()
                ch,plane,wire = [int(x) for x in chunks[:3]]
                beg = tuple([float(x) for x in chunks[3:6]])
                end = tuple([float(x) for x in chunks[6:9]])
                wire = Wire(ch, plane, wire, beg, end)
                ymax = max(beg[1], end[1])
                wires.append(((plane, ymax), wire))
        wires.sort()                          # by plane/ymax
        self.wires = [w[1] for w in wires]
        return self.wires

    @property
    def bounding_box(self):
        mins = list()
        maxs = list()

        for axis in range(3):
            pts = list()
            pts += [w.beg[axis] for w in self.wires]
            pts += [w.end[axis] for w in self.wires]
            mins.append(min(pts))
            maxs.append(max(pts))
        return (tuple(mins), tuple(maxs))
        

def plot_wires(wobj, wire_filter=None):
    bbmin, bbmax = wobj.bounding_box
    xmin, xmax = bbmin[2],bbmax[2]
    ymin, ymax = bbmin[1],bbmax[1]
    dx = xmax-xmin
    dy = ymax-ymin
    wires = wobj.wires

    print (xmin,ymin), (dx,dy)

    wirenums = [w.wire for w in wires]
    minwire = min(wirenums)
    maxwire = max(wirenums)
    nwires = maxwire-minwire+1

    if wire_filter:
        wires = [w for w in wires if wire_filter(w)]
    ax = plt.axes()
    ax.set_aspect('equal', 'box') #'datalim')
    ax.add_patch(mpatches.Rectangle((xmin, ymin), dx, dy,
                                    color="black", fill=False))

    cmap = plt.get_cmap('rainbow')        # seismic is bluewhitered

    colors = [cmap(i) for i in numpy.linspace(0, 1, nwires)]
    for ind, one in enumerate(wires):
        color = colors[one.wire-minwire]
        x = numpy.asarray((one.beg[2], one.end[2]))
        y = numpy.asarray((one.beg[1], one.end[1]))
        plt.plot(x, y, color=color)

    plt.axis([xmin,xmax,ymin,ymax])
    
