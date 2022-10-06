import simfish
import pandas as pd
import random
import numpy as np
from skimage.transform import resize
from skimage import io
import matplotlib.pyplot as plt
from tqdm import tqdm

def create_spot_image():
    imgs = simfish.simulate_images(6, 2, n_spots = 50, random_n_spots=True, image_shape=(100,100), sigma=15, noise_level=0) 

    fig, axs = plt.subplots(3,2)
    axs = axs.flatten()
    for i, (img, ground_truth) in enumerate(imgs):
        print(ground_truth.shape)
        axs[i].imshow(img)
    plt.show()

def resizePointWithTransform(point, zoom_factors):
     """resizePointWithTransform.                      
 
     Parameters
     ----------
     point :
         2D point in form of a tuple
     zoom_factors :
         transform in the form of the zoom factors used in the ndi.zoom function by scipy
     """  
     point = np.array([point]) # add a second artificial dimension to it
     new_point = np.array(zoom_factors) * point
     new_point = (float(new_point[0,0]), float(new_point[0,1]))
     return new_point


def createGeneExpressionImg(rows, cols, shape=(100,100)):
    zeros = np.zeros(shape, dtype=np.float)
    for row,col in zip(rows, cols):
        zeros[int( row ),int( col )] = 1.0
    return zeros

if __name__ == '__main__':
    # simfish.load_extract_template(".")
    pattern_strings = ["random" , "foci", "intranuclear", "extranuclear", "nuclear_edge", "perinuclear", "cell_edge", "pericellular", "protrusion"]
    for pattern in tqdm( pattern_strings ):
        for i in tqdm(range(317)):
            for it in range(40):
                n_spots = random.randint(600, 1100)
                try:
                    instance_coord = simfish.simulate_localization_pattern("/media/yob/simulated_cells/templates/", i_cell = i, n_spots = n_spots, pattern = pattern, proportion_pattern = 0.9)
                except ValueError:
                    continue
                rows = instance_coord["rna_coord"][:,1]
                cols = instance_coord["rna_coord"][:,2]
                zoom_factors = 100 / instance_coord["nuc_mask"].shape[0], 100 / instance_coord["nuc_mask"].shape[1]
                new_nuc = resize(instance_coord["nuc_mask"], (100,100))
                new_cell = resize(instance_coord["cell_mask"], (100,100))
                new_rows = []
                new_cols = []

                for point in zip(rows, cols):
                    new_point = resizePointWithTransform(point, zoom_factors)
                    new_rows.append(new_point[0])
                    new_cols.append(new_point[1])

                df = pd.DataFrame(list(zip(new_rows, new_cols)), columns=["row", "col"])


                ge_image = createGeneExpressionImg(new_rows, new_cols, (100,100))

                ## save:
                # io.imsave(f"/media/yob/simulated_cells/nuclei/{pattern}_{it}it_{n_spots}spots_icell{2}_nucleus.tif", new_nuc, check_contrast=False)
                # io.imsave(f"/media/yob/simulated_cells/cell/{pattern}_{it}it_{n_spots}spots_icell{2}_cell.tif", new_cell, check_contrast=False)
                io.imsave(f"/media/yob/simulated_cells/actually_all_cells_non-augmented/gene_expression/{pattern}_{it}it_{n_spots}spots_icell{i}_{it}it_gene_expression.tif", ge_image, check_contrast=False)
                df.to_csv(f"/media/yob/simulated_cells/actually_all_cells_non-augmented/gene_coords/{pattern}_{it}it_{n_spots}spots_icell{i}_{it}it_coords.tif", index=False)

