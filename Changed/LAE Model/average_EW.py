import numpy as np

EW = lambda x: 6.67298*x +80.96867
EWs = [EW(i) for i in np.linspace(0,14,10000)]
print(np.mean(EWs))