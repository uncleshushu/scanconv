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

imgs = dict.fromkeys(['bicubic', 'bicubic_ocl','bilinear_ocl', 'nearest_ocl'], None)

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


fig_2 = plt.figure('evaluation', figsize=(12, 4))
mse_plt = plt.subplot2grid((1, 3), (0, 0))
mse_plt.set_title('MSE')
mse_plt.bar(mses.keys(), mses.values(), alpha=0.8)
for k,v in mses.items():    
    mse_plt.text(k, v, '%.4f' % v, ha='center', va= 'bottom')

psnr_plt = plt.subplot2grid((1, 3), (0, 1))
psnr_plt.set_title('PSNR')
# psnr_plt.xlabel('method')
psnr_plt.bar(psnrs.keys(), psnrs.values(), alpha=0.8)
for k,v in psnrs.items():    
    psnr_plt.text(k, v, '%.2f' % v, ha='center', va= 'bottom')

ssim_plt = plt.subplot2grid((1, 3), (0, 2))
ssim_plt.set_title('SSIM')
ssim_plt.bar(ssims.keys(), ssims.values(), alpha=0.8)
for k,v in ssims.items():    
    ssim_plt.text(k, v, '%.2f' % v, ha='center', va= 'bottom')

# fig_1.show()
# fig_2.show()
plt.show()