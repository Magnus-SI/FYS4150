import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.figure()
for N in ('10', '100', '1000'):
    df = pd.read_csv('N' + N + '.csv')
    x = df['x'].values
    v = df['v'].values
    plt.plot(x, v, label = 'N = %s'%N)

u = df['u'].values
plt.plot(x, u, label = 'closed form sol.')
plt.legend()
plt.show()
