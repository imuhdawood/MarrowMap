import os
from openslide import OpenSlide
# import time
import argparse
import gc
import numpy as np
from shutil import rmtree
from tempfile import mkdtemp
import cv2
import math
import tifffile

######################################################
parser = argparse.ArgumentParser()

parser.add_argument('--image_folder', required=True, type=str, help='Path to the whole slide image')
parser.add_argument('--mask_folder', required=True, type=str, help='Path to the tissue mask')
parser.add_argument('--tile_size', default=512, type=int, help='The tile size')
parser.add_argument('--stride_ratio', default=4, type=int, help='The stride ratio')
parser.add_argument('--result_folder', required=True, type=str, help='Output path')

args = parser.parse_args()


def main():
    
    round_grp = ['Xenium']
    
    for r in range(len(round_grp)):
        
        sub_path = '{}/{}'.format(args.result_folder,round_grp[r])
        
        filelist=os.listdir('{}/{}'.format(args.image_folder,round_grp[r]))
        for i in range(len(filelist)):
            
            if os.path.exists('{}/{}_lesswhite.ome.tif'.format(sub_path,filelist[i].split('.ndpi')[0])):
                continue
            
            filename = '{}/{}/{}'.format(args.image_folder,round_grp[r],filelist[i])
            if filename.split('/')[-1].endswith('.ndpi'):
                imgname = filename.split('/')[-1].split('.ndpi')[0]
            else:
                imgname = filename.split('/')[-1].split('.ome.tif')[0]
            
            if not os.path.exists('{}/{}/progress_coord.txt'.format(sub_path,imgname)):
                print('{} does not exist'.format(imgname))
                continue
            else:
                if not os.path.exists('{}/{}/tile_scores.txt'.format(sub_path,imgname)):
                    print('Tile score for {} does not exist'.format(imgname))
                    continue
            
            print('Started {}'.format(imgname))
            
            wsi_img = OpenSlide(filename)
            
            sizes_xy = wsi_img.dimensions
            
            maskname_cell = '{}/{}/{}.ome.tif'.format(args.mask_folder,round_grp[r],imgname)
            if not os.path.exists(maskname_cell):
                print('{} does not exist'.format(filename.split('/')[-1].split('.')[0].split('_')[-1]))
                continue
            
            pixsize_ratio = round(1.0)
            
            with open('{}/{}/progress_coord.txt'.format(sub_path,imgname)) as p:
                strings_p = p.readlines()
                progress_coord  = [string.split(' ') for string in strings_p]
                p.close()
            
            tile_scores=[]
            
            with open('{}/{}/tile_scores.txt'.format(sub_path,imgname)) as p:
                strings_p = p.readlines()
                tile_scores  = [string.split(' ') for string in strings_p]

                p.close()
            
            coordinates=[]
            coordinates_x=[]
            coordinates_y=[]
            
            orig_idx = []
            
            for ii in range(len(progress_coord)-1):
                for iii in range(2,len(progress_coord[ii]),2):
                    coordinates.append((int(progress_coord[ii][iii+1]),int(progress_coord[ii][iii])))
                    orig_idx.append(ii)
                    coordinates_x.append(int(progress_coord[ii][iii]))
                    coordinates_y.append(int(progress_coord[ii][iii+1]))
                
            # Find and fix unique_coord
            # This is to avoid the case where the difference in coordinate is 1.
            # If there is the case, change the larger to smaller value
            # exmaple: (10,10) (11,20) --> (10,10) (10,20)
            
            coord_x_fix = np.copy(coordinates_x)
            coord_y_fix = np.copy(coordinates_y)
            
            sort_idx_x = np.argsort(coord_x_fix)
            sort_idx_y = np.argsort(coord_y_fix)
            
            sorted_x = coord_x_fix[sort_idx_x]
            sorted_y = coord_y_fix[sort_idx_y]
            
            diff_x = np.diff(sorted_x)
            diff_y = np.diff(sorted_y)
            
            one_diff_idx_x = np.where(diff_x==1)[0]
            one_diff_idx_y = np.where(diff_y==1)[0]
            
            
            for iii in one_diff_idx_x:
                cnt_idx = 0
                base_v = sorted_x[iii]
                while diff_x[iii+1+cnt_idx] == 0:
                    coord_x_fix[sort_idx_x[iii+1+cnt_idx]] = base_v
                    cnt_idx+=1
                coord_x_fix[sort_idx_x[iii+1+cnt_idx]] = base_v
            
            for iii in one_diff_idx_y:
                cnt_idx = 0
                base_v = sorted_y[iii]
                while diff_y[iii+1+cnt_idx] == 0:
                    coord_y_fix[sort_idx_y[iii+1+cnt_idx]] = base_v
                    cnt_idx+=1
                coord_y_fix[sort_idx_y[iii+1+cnt_idx]] = base_v
            
            # pdb.set_trace()
            
            del sort_idx_x
            del sort_idx_y
            del sorted_x
            del sorted_y
            del diff_x
            del diff_y
            del one_diff_idx_x
            del one_diff_idx_y
            
            coord_fix = np.copy(coordinates)
            
            for iii in range(len(coord_fix)):
                coord_fix[iii][1] = coord_x_fix[iii]
                coord_fix[iii][0] = coord_y_fix[iii]
            
            unique_coord = np.unique(coord_fix, axis=0)
            
            temp_folder1 = mkdtemp(dir='../code')
            temp_filename1 = os.path.join(temp_folder1, 'file.dat')
                
            canvas1 = np.memmap(temp_filename1,
                                dtype='float32',
                                mode='w+',
                                shape=(int(sizes_xy[1]/10), int(sizes_xy[0]/10)))            
            
            overlap_score_txt=open('{}/{}/overlap_scores.txt'.format(sub_path,imgname),'w')
            
            tissue_mask = tifffile.imread(maskname_cell)
            if (tissue_mask.shape[1] != int(sizes_xy[0]/10)) | (tissue_mask.shape[0] != int(sizes_xy[1]/10)):
                tissue_mask = cv2.resize(tissue_mask,(int(sizes_xy[0]/10),int(sizes_xy[1]/10)),interpolation=cv2.INTER_NEAREST)
            
            for ii in range(len(unique_coord)):
                idx = np.where((np.array(coord_x_fix)==unique_coord[ii][1])&(np.array(coord_y_fix)==unique_coord[ii][0]))[0]
                
                if len(idx) > 1:
                    score_gather = []
                    for iii in range(len(idx)):
                        if tile_scores[orig_idx[idx[iii]]][0] != 'None\n':
                            score_gather.append(float(tile_scores[orig_idx[idx[iii]]][0]))
                    if len(score_gather) == 0:
                        final_tile_score = np.nan
                    else:
                        final_tile_score = np.mean(score_gather)
                else:
                    if tile_scores[orig_idx[idx[0]]][0] != 'None\n':
                        final_tile_score = float(tile_scores[orig_idx[idx[0]]][0])
                    else:
                        final_tile_score = np.nan
    
                if math.isnan(final_tile_score):
                    print('{}: {} out of {} has been done'.format(filelist[i],ii+1,len(unique_coord)))
                    continue
    
                strt_tile_x = unique_coord[ii][1]-int((args.tile_size/(args.stride_ratio*2))/pixsize_ratio)
                if strt_tile_x < 0:
                    strt_tile_x = 0
    
                strt_tile_y = unique_coord[ii][0]-int((args.tile_size/(args.stride_ratio*2))/pixsize_ratio)
                if strt_tile_y < 0:
                    strt_tile_y = 0
    
                end_tile_x = unique_coord[ii][1]+int((args.tile_size/(args.stride_ratio*2))/pixsize_ratio)
                if end_tile_x >= sizes_xy[0]:
                    end_tile_x = sizes_xy[0]-1
    
                end_tile_y = unique_coord[ii][0]+int((args.tile_size/(args.stride_ratio*2))/pixsize_ratio)
                if end_tile_y >= sizes_xy[1]:
                    end_tile_y = sizes_xy[1]-1
                
                if end_tile_x <= strt_tile_x:
                    continue
                
                if end_tile_y <= strt_tile_y:
                    continue
            
                crop_mask = tissue_mask[int(strt_tile_y/10):int(end_tile_y/10), int(strt_tile_x/10):int(end_tile_x/10)]
                if (((crop_mask == 0).sum() / crop_mask.size) >= 0.75):
                    print('{} out of {} has been done'.format(ii+1,len(unique_coord)))
                    continue
                
                del crop_mask
                
                overlap_score_txt.write('X:{} Y:{} Score:{}\n'.format(strt_tile_x,strt_tile_y,final_tile_score))
                canvas1[int(strt_tile_y/10):int(end_tile_y/10)+1, int(strt_tile_x/10):int(end_tile_x/10)+1] = final_tile_score
                print('{} out of {} has been done'.format(ii+1,len(unique_coord)))
                
            overlap_score_txt.write('finished')
            overlap_score_txt.close()
            savename = '{}/{}_lesswhite.ome.tif'.format(sub_path, imgname)
            cv2.imwrite(savename, canvas1)
            print('{} has been done'.format(imgname))
            del canvas1
            
            gc.collect()
            rmtree(temp_folder1)
            
if __name__ == '__main__':
    main()
