import json
import matplotlib.pyplot as plt
import numpy as np

# loading json file
with open('conf.json', 'r') as file:
    data = json.load(file)


pae_data = data.get('pae', [])


pae_matrix = np.array(pae_data)


plt.style.use('seaborn-whitegrid')


plt.figure(figsize=(10, 8))  # adjust size
plt.imshow(pae_matrix, cmap='viridis', interpolation='nearest')  

# add colorbar
plt.colorbar(label='PAE Value')

# add title
plt.title('Pairwise Affinity Error Grid Heatmap')

# x-y label
plt.xlabel('Residue Index')
plt.ylabel('Residue Index')

# x-y ticks
x_ticks = np.linspace(0, pae_matrix.shape[1] - 1, 5).astype(int)  
y_ticks = np.linspace(0, pae_matrix.shape[0] - 1, 5).astype(int)  


plt.xticks(x_ticks)  
plt.yticks(y_ticks)  
plt.grid(False)


# save as PDF file
plt.savefig('conf5.pdf', format='pdf')
