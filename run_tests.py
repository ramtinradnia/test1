from gasnetwork.typecurve.synthetictypecurve import SyntheticTypeCurve
from gasnetwork.utils import engunits
import matplotlib.pyplot as plt
import numpy as np
import array as arr
import json

#C:\Users\username\storagelocation\gasnetworkfinal\gasnetwork-master>python ./run_tests.py

q_zero = 0 # MMSCFD
q_zero = engunits.convert_value('svf', q_zero, 'MMSCFD')

q_peak = 10 # MMSCFD
q_peak = engunits.convert_value('svf', q_peak, 'MMSCFD')

type_curve = SyntheticTypeCurve(q_zero, q_peak, 5, 6, 0, 1, 1)
ts = list(range(12 * 30))
qs = type_curve.build_curve(ts)
qs_conv = [None]*len(qs)

for i in range(len(qs)):
    qs_conv[i] = engunits.convert_value('svf', qs[i], 'MMSCFD')
print(qs_conv)

print (ts)  # in months
print (qs)  # in MMSCFD

fig = plt.plot(ts, qs, 'b')

ref_series = json.load(open('ref_curve.py'))
qs=ref_series["qs"]

    
ts=ref_series["ts"]
print(ts)
print(qs)
plt.plot(ts, qs, 'r')

plt.show()




