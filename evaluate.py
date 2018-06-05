#! /usr/bin/env python3

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from skimage.measure import compare_mse, compare_psnr, compare_ssim
import math

def mse(img1, img2):
    return np.mean( (img1 - img2) ** 2 )

def psnr(img1, img2):
    mse = np.mean( (img1 - img2) ** 2 )
    if mse == 0:
        return 100
    PIXEL_MAX = 255.0
    return 20 * math.log10(PIXEL_MAX / math.sqrt(mse))

imgs = dict.fromkeys(['bicubic', 'bicubic_ocl', 'col_linear_row_cubic_ocl', 'bilinear_ocl', 'nearest_ocl'], None)

for k in imgs.keys():
    imgs[k] = mpimg.imread(k+'.png')


mses = {}
psnrs = {}
ssims = {}
for k,v in imgs.items():
    if k != 'bicubic':
        mses[k] = compare_mse(imgs['bicubic'], v)
        psnrs[k] = compare_psnr(imgs['bicubic'], v)
        ssims[k] = compare_ssim(imgs['bicubic'], v)

print("MSE:", mses)
print("PSNR:", psnrs)
print("SSIM", ssims)

fig_1 = plt.figure('results', figsize=(15, 9))

row_cnt = math.ceil(len(imgs)/3.0)
for i, (k, v) in enumerate(imgs.items()):
    p = plt.subplot2grid((row_cnt, 3), (i//3, i%3))
    p.set_title(k)
    p.imshow(v,cmap='Greys_r')
    p.axis('off')

fig_1.tight_layout()

fig_2 = plt.figure('evaluation', figsize=(9, 9))

measures = {'MSE': mses, 'PSNR': psnrs, 'SSIM': ssims}
for i, (name, data) in enumerate(measures.items()):
    ax = plt.subplot2grid((3, 1), (i, 0))
    ax.set_title(name)
    ax.barh(list(data.keys()), list(data.values()), alpha=0.8)
    # ax.barh(range(len(data)), list(data.values()), alpha=0.8, tick_label=data.keys())
    for k,v in data.items():    
        ax.text(v, k, '%.3g' % v, ha='left', va= 'center')

fig_2.tight_layout()

# fig_1.show()
# fig_2.show()
plt.show()