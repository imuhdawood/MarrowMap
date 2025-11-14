import os
from openslide import OpenSlide
import argparse
import gc
import numpy as np
import model
import torch
import torch.backends.cudnn as cudnn
import torchvision
import tqdm
import tifffile
import cv2

######################################################
parser = argparse.ArgumentParser()

parser.add_argument('--image_folder', required=True, type=str, help='Path to the whole slide image')
parser.add_argument('--mask_folder', required=True, type=str, help='Path to the tissue mask')
parser.add_argument('--result_folder', required=True, type=str, help='Output path')
parser.add_argument('--resume_fib', required=True, type=str, help='Path to the CIF model')
parser.add_argument('--tile_magnification', default=40, type=int, help='The magnification which the tile needs to be extracted')
parser.add_argument('--wsi_highest_magnification', default=40, type=int, help='Set the magnification of whole slide image')
parser.add_argument('--tile_size', default=512, type=int, help='The tile size')
parser.add_argument('--stride_ratio', default=4, type=int, help='The stride')
parser.add_argument('--norm', default='y', type=str, help='Whether to use normalization or not')
parser.add_argument('--batch_size', default=1, type=int, help='The batch size')
parser.add_argument('--workers', default=0, type=int, help='The number of workers')
parser.add_argument('--gpu_id', default='0', type=str, help='To specfify the GPU IDs')

args = parser.parse_args()

class TileDataset(torch.utils.data.Dataset): # Get the tile set from whole slide image

    def __init__(self, reader, mask_cell_reader, img_dim, 
                 roi_coords, img_name, wsi_highest_magnification, tile_magnification, tile_size, txt_resultpath, 
                 pix_ratio, stride_ratio, norm):
        self.reader = reader
        self.mask_cell_reader = mask_cell_reader
        self.img_dim = img_dim
        self.roi_coords = roi_coords
        self.img_name = img_name
        self.wsi_highest_magnification = wsi_highest_magnification
        self.tile_magnification = tile_magnification
        self.tile_size = tile_size
        self.txt_resultpath = txt_resultpath
        self.stride_ratio = stride_ratio
        self.pix_ratio = pix_ratio
        self.norm = norm
        self.stride = int(self.tile_size*(self.wsi_highest_magnification / self.tile_magnification)/self.stride_ratio)
        self.rescale_tile_size, 
        self.rescale_factor, 
        self.downsample = self._get_level()
        self.coordinates_list= self._get_coordinates()
        self.totensor = torchvision.transforms.ToTensor()
        self.normalise = torchvision.transforms.Normalize((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))

    def _get_level(self):

        downsample = self.wsi_highest_magnification / self.tile_magnification
        rescale_tile_size = int(np.round(self.tile_size * downsample))
        rescale_factor = 1/downsample
        
        return rescale_tile_size, rescale_factor, downsample

    def _get_coordinates(self): # Find end pixel coordinates

        coordinates = []
        
        if os.path.exists('{}/{}/progress_coord.txt'.format(self.txt_resultpath,self.img_name)):
            with open('{}/{}/progress_coord.txt'.format(self.txt_resultpath,self.img_name)) as p:
                strings_p = p.readlines()
                progress_coord  = [string.split(' ') for string in strings_p]
                p.close()
            
            if progress_coord[-1][0] == 'finished':
                for ii in range(len(progress_coord)-1):
                    coordinates.append((int(progress_coord[ii][1].split(':')[1]),int(progress_coord[ii][0].split(':')[1])))
                print('Coord of {} already exists'.format(self.img_name))
                return coordinates
        
        start_row = self.roi_coords[2]
        start_col = self.roi_coords[0]
        end_row = self.roi_coords[3]
        end_col = self.roi_coords[1]
        
        coordinates = []
        
        if not os.path.exists('{}/{}'.format(self.txt_resultpath,self.img_name)):
            os.makedirs('{}/{}'.format(self.txt_resultpath,self.img_name))
        
        if os.path.exists('{}/{}/progress_coord.txt'.format(self.txt_resultpath,self.img_name)):
            with open('{}/{}/progress_coord.txt'.format(self.txt_resultpath,self.img_name)) as p:
                strings_p = p.readlines()
                progress_coord  = [string.split(' ') for string in strings_p]
                p.close()
            
            if len(progress_coord[-1]) < 3:
                for ii in range(len(progress_coord)-1):
                    coordinates.append((int(progress_coord[ii][1].split(':')[1]),int(progress_coord[ii][0].split(':')[1])))
                start_row = int(progress_coord[-2][1].split(':')[-1]) + int(self.tile_size/self.pix_ratio)
                if start_row > self.img_dim[1]:
                    start_col = int(progress_coord[-2][0].split(':')[-1]) + int(self.tile_size/self.pix_ratio)
                else:
                    start_col = int(progress_coord[-2][0].split(':')[-1])
            else:
                for ii in range(len(progress_coord)):
                    coordinates.append((int(progress_coord[ii][1].split(':')[1]),int(progress_coord[ii][0].split(':')[1])))
                start_row = int(progress_coord[-1][1].split(':')[-1]) + int(self.tile_size/self.pix_ratio)
                if start_row > self.img_dim[1]:
                    start_col = int(progress_coord[-1][0].split(':')[-1]) + int(self.tile_size/self.pix_ratio)
                    start_row = 0
                else:
                    start_col = int(progress_coord[-1][0].split(':')[-1])
            
            print('Continue getting the coordinates from X:{} Y:{}'.format(start_col,start_row))
        
        print('Start getting the coordinates')
        
        cnt_trial = 0
        
        len_trial_col = len(np.arange(start_col, end_col, int((self.tile_size/self.pix_ratio)/self.stride_ratio)))
        len_trial_row = len(np.arange(start_row, end_row, int((self.tile_size/self.pix_ratio)/self.stride_ratio)))
        
        for col in np.arange(start_col, end_col, int((self.tile_size/self.pix_ratio)/self.stride_ratio)):
            for row in np.arange(start_row, end_row, int((self.tile_size/self.pix_ratio)/self.stride_ratio)):
                X = int(col)
                Y = int(row)
                W = int(self.tile_size/self.pix_ratio*self.downsample)
                H = int(self.tile_size/self.pix_ratio*self.downsample)
                if X + W > self.img_dim[0]:
                    W = self.img_dim[0] - X
                    if W < int(341/self.pix_ratio):
                        cnt_trial+=1
                        continue
                if Y + H > self.img_dim[1]:
                    H = self.img_dim[1] - Y
                    if H < int(341/self.pix_ratio):
                        cnt_trial+=1
                        continue
                
                cnt_trial+=1
                
                if cnt_trial >= len_trial_col*len_trial_row:
                    if not os.path.exists('{}/{}/progress_coord.txt'.format(self.txt_resultpath,self.img_name)):
                        txt_file = open('{}/{}/progress_coord.txt'.format(self.txt_resultpath,self.img_name), 'x')
                    else:
                        txt_file = open('{}/{}/progress_coord.txt'.format(self.txt_resultpath,self.img_name), 'a')
                    txt_file.write('finished')
                    txt_file.close()
                
                
                crop_mask = self.mask_cell_reader[int(Y/10):int((Y+H)/10),int(X/10):int((X+W)/10)]
                
                if not os.path.exists('{}/{}/total_numb_tile.txt'.format(self.txt_resultpath,self.img_name)):
                    total_tile_txt = open('{}/{}/total_numb_tile.txt'.format(self.txt_resultpath,self.img_name), 'x')
                else:
                    total_tile_txt = open('{}/{}/total_numb_tile.txt'.format(self.txt_resultpath,self.img_name), 'a')
                
                total_tile_txt.write('X:{} Y:{} '.format(X, Y))
                
                for cnt_tile_x in range(1,int(self.stride_ratio*2),2):
                    for cnt_tile_y in range(1,int(self.stride_ratio*2),2):
                        if (cnt_tile_x == int(self.stride_ratio*2)-1) & (cnt_tile_y == int(self.stride_ratio*2)-1):
                            total_tile_txt.write('{} {}\n'.format(X+int((self.tile_size/(self.stride_ratio*2))/self.pix_ratio)*cnt_tile_x,
                                                                  Y+int((self.tile_size/(self.stride_ratio*2))/self.pix_ratio)*cnt_tile_y))
                        else:
                            total_tile_txt.write('{} {} '.format(X+int((self.tile_size/(self.stride_ratio*2))/self.pix_ratio)*cnt_tile_x,
                                                                 Y+int((self.tile_size/(self.stride_ratio*2))/self.pix_ratio)*cnt_tile_y))
                
                if (((crop_mask < 1).sum() / crop_mask.size) >= 0.5):
                    continue
                
                print('X: {} Y: {}'.format(X,Y))
                if not os.path.exists('{}/{}/progress_coord.txt'.format(self.txt_resultpath,self.img_name)):
                    txt_file = open('{}/{}/progress_coord.txt'.format(self.txt_resultpath,self.img_name), 'x')
                else:
                    txt_file = open('{}/{}/progress_coord.txt'.format(self.txt_resultpath,self.img_name), 'a')
                
                txt_file.write('X:{} Y:{} '.format(X, Y))
                
                for cnt_tile_x in range(1,int(self.stride_ratio*2),2):
                    for cnt_tile_y in range(1,int(self.stride_ratio*2),2):
                        if (cnt_tile_x == int(self.stride_ratio*2)-1) & (cnt_tile_y == int(self.stride_ratio*2)-1):
                            txt_file.write('{} {}\n'.format(X+int((self.tile_size/(self.stride_ratio*2))/self.pix_ratio)*cnt_tile_x,
                                                            Y+int((self.tile_size/(self.stride_ratio*2))/self.pix_ratio)*cnt_tile_y))
                        else:
                            txt_file.write('{} {} '.format(X+int((self.tile_size/(self.stride_ratio*2))/self.pix_ratio)*cnt_tile_x,
                                                           Y+int((self.tile_size/(self.stride_ratio*2))/self.pix_ratio)*cnt_tile_y))
                
                coordinates.append((Y,X))
                
                
                txt_file.close()
        
        total_tile_txt.close()
        print('Finished getting the coordinates')
        return coordinates

    def __len__(self):
        return len(self.coordinates_list)

    def __getitem__(self, idx):
        
        row, col = self.coordinates_list[idx]
        X = col
        Y = row
        W = int(self.rescale_tile_size/self.pix_ratio)
        H = int(self.rescale_tile_size/self.pix_ratio)
        
        if (X+W > self.img_dim[0]) & (Y+H > self.img_dim[1]):
            W_temp = self.img_dim[0]-X
            H_temp = self.img_dim[1]-Y
            image = self.reader.read_region((col, row), 0, (W_temp,H_temp))
            end_x = self.img_dim[0]
            end_y = self.img_dim[1]
        elif (X+W > self.img_dim[0]) & (Y+H < self.img_dim[1]):
            W_temp = self.img_dim[0]-X
            image = self.reader.read_region((col, row), 0, (W_temp,int(self.rescale_tile_size/self.pix_ratio)))
            end_x = self.img_dim[0]
            end_y = 0
        elif (X+W < self.img_dim[0]) & (Y+H > self.img_dim[1]):
            H_temp = self.img_dim[1]-Y
            image = self.reader.read_region((col, row), 0, (int(self.rescale_tile_size/self.pix_ratio), H_temp))
            end_x = 0
            end_y = self.img_dim[1]
        else:
            image = self.reader.read_region((col, row), 0, (int(self.rescale_tile_size/self.pix_ratio),int(self.rescale_tile_size/self.pix_ratio)))

            end_x = 0
            end_y = 0
            end_x = 0
            end_y = 0
        
        image = np.array(image.convert("RGB"))
        
        image_orig = cv2.resize(image, (512, 512), interpolation=cv2.INTER_AREA)
        
        if self.norm == 'y':
            
            clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8))
            HSV = cv2.cvtColor(image_orig, cv2.COLOR_RGB2HSV)
            HSV[:, :, 0] = clahe.apply(HSV[:, :, 0])
            image = cv2.cvtColor(HSV, cv2.COLOR_HSV2RGB)
            
            image = image.astype(np.uint8)
        
        else:
            image = image_orig.astype(np.uint8)

        # convert to tensor and normalise
        image = self.normalise(self.totensor(image))
        
        # return an image and the corresponding coordinates
        return image, (row, col), (end_y,end_x)

    def _collate_filter_out_none(self, batch):
        "Puts each data field into a tensor with outer dimension batch size"
        batch = [x for x in batch if x is not None]
        if len(batch):
            return torch.utils.data.dataloader.default_collate(batch)
        else:
            return None, (None, None)



def main():
    # initialize CUDA
    os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"  # see issue #152
    os.environ["CUDA_VISIBLE_DEVICES"] = args.gpu_id
    cudnn.benchmark = True
    
    os.environ['KMP_DUPLICATE_LIB_OK']='True'
    
    #------------------------------------------------------------------------
    # create model
    net = model.Siamese_rank_bn_noGAP_HE(channel_choice = 3)
    net = torch.nn.DataParallel(net)

    if torch.cuda.is_available():
        net = net.cuda()
    
    if args.resume_fib:
        if os.path.isfile(args.resume_fib):
            print("=> loading checkpoint '{}'".format(args.resume_fib))
            checkpoint = torch.load(args.resume_fib)
            net.load_state_dict(checkpoint['state_dict'])

            print("=> loaded checkpoint '{}'"
                  .format(args.resume_fib))
        else:
            raise IOError("=> no checkpoint found at '{}'".format(args.resume_fib))
    
    #------------------------------------------------------------------------
    
    round_grp = ['Xenium'] # Specify the project/group name
    
    # The configuration of normalization of the CIF score -------------------
    denom = 45.9749145507812
    
    numer = -19.0318393707275
    
    #------------------------------------------------------------------------
    
    for r in range(len(round_grp)):
        sub_path = '{}/{}'.format(args.result_folder,round_grp[r])
        
        if not os.path.exists(sub_path):
            os.makedirs(sub_path)
            
        filelist=os.listdir('{}/{}'.format(args.image_folder,round_grp[r]))
        for i in range(len(filelist)):
            
            filename = '{}/{}/{}'.format(args.image_folder,round_grp[r],filelist[i])
            if filename.endswith('.ndpi'):
                imgname = filename.split('/')[-1].split('.ndpi')[0]
            else:
                imgname = filename.split('/')[-1].split('.ome')[0]
            
            if os.path.exists('{}/{}/progress_coord.txt'.format(sub_path,imgname)):
                with open('{}/{}/progress_coord.txt'.format(sub_path,imgname)) as p:
                    strings_p = p.readlines()
                    progress_coord  = [string.split(' ') for string in strings_p]
                    p.close()
                
                if progress_coord[-1][0] == 'finished':
                    if os.path.exists('{}/{}/tile_scores.txt'.format(sub_path,imgname)):
                        print('{} already exists'.format(imgname))
                        continue
            
            maskname_cell = '{}/{}/{}.ome.tif'.format(args.mask_folder,round_grp[r],imgname)
            
            if not os.path.exists(maskname_cell):
                print('{} does not exist'.format(imgname))
                continue
            # pdb.set_trace()
            
            print('Started {}'.format(filename))
            
            wsi_img = OpenSlide(filename)
            
            sizes_xy = wsi_img.dimensions
            
            pixsize_ratio = round(1.0)
            
            tissue_mask = tifffile.imread(maskname_cell)
            
            if (tissue_mask.shape[1] != int(sizes_xy[0]/10)) | (tissue_mask.shape[0] != int(sizes_xy[1]/10)):
                tissue_mask = cv2.resize(tissue_mask,(int(sizes_xy[0]/10),int(sizes_xy[1]/10)),interpolation=cv2.INTER_NEAREST)
            
            tissue_idx = np.where(np.array(tissue_mask) > 0)
            
            strt_x = np.min(tissue_idx[1])*10
            
            strt_y = np.min(tissue_idx[0])*10
            
            end_x = np.max(tissue_idx[1])*10
            
            end_y = np.max(tissue_idx[0])*10
            
            if end_x >= sizes_xy[0]:
                end_x = sizes_xy[0]-1
            
            if end_y >= sizes_xy[1]:
                end_y = sizes_xy[1]-1
            
            if strt_x < 0:
                strt_x = 0
            
            if strt_y < 0:
                strt_y = 0
            
            roi_region = [strt_x, end_x, strt_y, end_y]
             
            dataset = TileDataset(reader = wsi_img,
                                 mask_cell_reader = tissue_mask,
                                 img_dim = sizes_xy,
                                 roi_coords = roi_region,
                                 img_name = imgname,
                                 wsi_highest_magnification=args.wsi_highest_magnification,
                                 tile_magnification=args.tile_magnification,
                                 tile_size=args.tile_size,
                                 txt_resultpath = sub_path,
                                 pix_ratio = pixsize_ratio,
                                 stride_ratio = args.stride_ratio,
                                 norm = args.norm)
            
            data_loader = torch.utils.data.DataLoader(dataset,
                                                      batch_size=args.batch_size,
                                                      shuffle=False,
                                                      num_workers=args.workers,
                                                      collate_fn=dataset._collate_filter_out_none)
            
            gc.collect()

            score_txt=open('{}/{}/tile_scores.txt'.format(sub_path,imgname),'w')
 
            # tile inference
            for images, (rows, cols), (end_ys, end_xs) in tqdm.tqdm(data_loader):
                if images is None:
                    continue

                if torch.cuda.is_available():
                    images = images.cuda()

                net.eval()

                with torch.no_grad():
                    pred_score = net(images)
                    
                pred_score_norm = (pred_score-numer)/denom
                pred_score_norm = pred_score_norm.cpu().numpy()
            
                if pred_score_norm > 1:
                    pred_score_norm = np.array([[1]])
                    pred_score_norm.astype(np.float32)
    
                elif pred_score_norm < 0:
                    pred_score_norm = np.array([[0]])
                    pred_score_norm.astype(np.float32)
    
                score_txt.write('{}\n'.format(pred_score_norm[0][0]))
    
                rows = rows.cpu().data.numpy()
                cols = cols.cpu().data.numpy()

                end_xs = end_xs.cpu().data.numpy()
                end_ys = end_ys.cpu().data.numpy()

                gc.collect()
            
            score_txt.close()
            gc.collect()
            

if __name__ == '__main__':
    main()
