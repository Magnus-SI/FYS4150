import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.figure()
for N in ('1', '2', '3'):
    df = pd.read_csv('values' + N + '.csv')
    x = df['x'].values
    v = df['v'].values
    plt.plot(x, v, label = 'N = %s'%N)

u = df['u'].values
plt.plot(x, u, label = 'closed form sol.')
plt.legend()
plt.show()
