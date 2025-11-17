import numpy as np
import matplotlib.pyplot as plt
import cv2
import os
import tifffile

img_path = 'path/to/data'
mask_path = 'path/to/CIFmask'
result_path = 'path/to/output'

if not os.path.exists(result_path):
    os.makedirs(result_path)


round_grp = ['Project Name']


for r in range(len(round_grp)):
    sub_img_path = '{}/downsampled/{}'.format(img_path,round_grp[r])
    sub_mask_path = '{}/{}'.format(mask_path,round_grp[r])
    sub_savepath = '{}/{}/heatmap'.format(result_path,round_grp[r])
    
    if not os.path.exists(sub_savepath):
        os.makedirs(sub_savepath)
    
    filelist = []
    
    for file in os.listdir(sub_mask_path):
        if file.endswith(".tif"):
            filelist.append(file.split('.ome.tif')[0])
        if file.endswith(".ndpi"):
            filelist.append(file.split('.ndpi')[0])
    
    
    for i in range(len(filelist)):
        imgname = filelist[i].split('_lesswhite')[0]
        
        imgid = imgname
        
        
        if os.path.exists('{}/{}.jpeg'.format(sub_savepath,imgname)):
            print('{} already exists'.format(imgname))
            continue

        aa = tifffile.imread('{}/{}_downsampled_x10.tif'.format(sub_img_path,imgname))
        if aa.shape[0] < 4:
            aa = np.moveaxis(aa,0,-1)
    
        img = aa
        img = img.astype(np.uint8)
        del aa
        
        mask = cv2.imread('{}/{}_lesswhite.ome.tif'.format(sub_mask_path,imgname),-1)
        
        img = cv2.resize(img,(mask.shape[1],mask.shape[0]))
        
        R, G, B = cv2.split(img)
        
        cmap = plt.get_cmap('jet')
        
        flag = mask > 0
        
        heatmap = (cmap(mask) * 255).astype(np.uint8)
        
        r, g, b = cv2.split(heatmap[:,:,0:3])

        R[flag] = (R[flag]*0.5) + (r[flag]*0.5)
        G[flag] = (G[flag]*0.5) + (g[flag]*0.5)
        B[flag] = (B[flag]*0.5) + (b[flag]*0.5)

        R = R.astype(np.uint8)
        G = G.astype(np.uint8)
        B = B.astype(np.uint8)

        overlay = cv2.merge([R,G,B])

        cv2.imwrite('{}/{}.jpeg'.format(sub_savepath,imgname),cv2.cvtColor(overlay, cv2.COLOR_BGR2RGB))
        print('{} has been finished'.format(imgname))
