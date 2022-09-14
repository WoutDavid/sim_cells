import matplotlib.pyplot as plt
from skimage import io
from glob import glob


img_list =glob( "/home/david/.config/nnn/mounts/nacho@10.38.76.144/yob/simulated_cells/gene_expression/*_1it_*spots_icell2_cell.tif")
pattern_strings = ["random" , "foci", "intranuclear", "extranuclear", "nuclear_edge", "perinuclear", "cell_edge", "pericellular", "protrusion"]


fig, axs = plt.subplots(2,5)
axs = axs.flatten()
for i, pattern in enumerate( pattern_strings ):
    pattenr_imgs = [file for file in img_list if pattern in file ]

    axs[i].imshow(io.imread(pattenr_imgs[0]))
    axs[i].set_title(pattern)
fig.delaxes(axs[-1])
plt.show()
