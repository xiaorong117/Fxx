#!/usr/bin/env python3
import math
start=42.73339209431653;end=start+5*0.01;time=start
for _ in range(5):time+=0.01
eps=32*math.ulp(max(1.0,abs(end)))
assert time<end and end-time<eps and not (time+eps<end)
print('PASS end-time floating tail is not treated as a physical timestep')
