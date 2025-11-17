import torch

# Convolutional layer with the option of using Instance-Batch Normalization (use_ibn = True)
# (Xingang Pan, Ping Luo, Jianping Shi, Xiaoou Tang. 
# "Two at Once: Enhancing Learning and Generalization Capacities via IBN-Net", ECCV2018)
# Note: This can also be used to train the model if the user wants to use either
#       Instance Normalization (use_instn = True) or Batch Normalization (use_bn = True).

class IBN(torch.nn.Module):
    def __init__(self, planes, ratio=0.5):
        super(IBN, self).__init__()
        self.half = int(planes * (1-ratio))
        self.ibn_BN = torch.nn.BatchNorm2d(self.half)
        self.ibn_IN = torch.nn.InstanceNorm2d(planes - self.half, affine=True)

    def forward(self, x):
        split = torch.split(x, self.half, 1)
        out1 = self.ibn_BN(split[0].contiguous())
        out2 = self.ibn_IN(split[1].contiguous())
        out = torch.cat((out1, out2), 1)
        return out

class convLayer_ibn(torch.nn.Module):
    def __init__(self, input_channel, output_channel, filter_size, stride, padding, 
                 pool_win, block, down_size, use_bn = False, use_instn = False, use_ibn = False):
        super(convLayer_ibn, self).__init__()
        if (use_bn == True) | (use_ibn == True):
            self.conv = torch.nn.Conv2d(input_channel, output_channel, filter_size, 
                                        stride=stride, padding = padding, bias = False)
            self.conv_same = torch.nn.Conv2d(output_channel, output_channel, filter_size, 
                                             stride=stride, padding = padding, bias = False)
            self.conv_same_1 = torch.nn.Conv2d(output_channel, output_channel, filter_size, 
                                               stride=stride, padding = padding, bias = False)
        else:
            self.conv = torch.nn.Conv2d(input_channel, output_channel, filter_size, 
                                        stride=stride, padding = padding)
            self.conv_same = torch.nn.Conv2d(output_channel, output_channel, filter_size, 
                                             stride=stride, padding = padding)
            self.conv_same_1 = torch.nn.Conv2d(output_channel, output_channel, filter_size, 
                                               stride=stride, padding = padding)
        
        self.use_bn = use_bn
        self.use_instn = use_instn
        self.use_ibn = use_ibn
        self.instn1 = torch.nn.InstanceNorm2d(output_channel)
        self.instn2 = torch.nn.InstanceNorm2d(output_channel)
        self.instn3 = torch.nn.InstanceNorm2d(output_channel)
        self.bn1 = torch.nn.BatchNorm2d(output_channel)
        self.bn2 = torch.nn.BatchNorm2d(output_channel)
        self.bn3 = torch.nn.BatchNorm2d(output_channel)
        self.IBN = IBN(output_channel, 0.5)
        self.max_pool = torch.nn.MaxPool2d(pool_win)
        self.relu = torch.nn.ReLU()
        self.down_size = down_size
        self.block = block
    
    def forward(self, x):
        if self.block > 1:
            for b in range(self.block):
                if b == 1:
                    if self.use_bn:
                        x = self.relu(self.bn2(self.conv_same(x)))
                    elif self.use_instn:
                        x = self.relu(self.instn2(self.conv_same(x)))
                    
                elif b == 2:
                    if self.use_bn:
                        x = self.relu(self.bn3(self.conv_same_1(x)))
                    elif self.use_instn:
                        x = self.relu(self.instn3(self.conv_same_1(x)))
                else:
                    if self.use_ibn:
                        x = self.relu(self.IBN(self.conv(x)))
                    elif self.use_bn:
                        x = self.relu(self.bn1(self.conv(x)))
                    elif self.use_instn:
                        x = self.relu(self.instn1(self.conv(x)))
        else:
            if self.use_bn:
                x = self.relu(self.bn1(self.conv(x)))
            elif self.use_instn:
                x = self.relu(self.instn1(self.conv(x)))
        if self.down_size:
            x = self.max_pool(x)
        return x


class Siamese_rank_bn_noGAP_HE(torch.nn.Module):
    def __init__(self, channel_choice):
        super(Siamese_rank_bn_noGAP_HE, self).__init__()
        
        if channel_choice == 3:
            self.down_block1 = convLayer_ibn(input_channel=3, output_channel=64, filter_size=3, 
                                             stride=1, padding=1, pool_win=(2,2), block=1, down_size=True,
                                             use_bn = True, use_instn = False, use_ibn = False)
        else:
            self.down_block1 = convLayer_ibn(input_channel=1, output_channel=64, filter_size=3, 
                                             stride=1, padding=1, pool_win=(2,2), block=1, down_size=True,
                                             use_bn = True, use_instn = False, use_ibn = False)
        
        self.down_block2 = convLayer_ibn(input_channel=64, output_channel=128, filter_size=3, 
                                     stride=1, padding=1, pool_win=(2,2), block=1, down_size=True,
                                     use_bn = True, use_instn = False, use_ibn = False)
        self.down_block3 = convLayer_ibn(input_channel=128, output_channel=256, filter_size=3, 
                                     stride=1, padding=1, pool_win=(2,2), block=1, down_size=True,
                                     use_bn = True, use_instn = False, use_ibn = False)
        self.down_block4 = convLayer_ibn(input_channel=256, output_channel=512, filter_size=3, 
                                     stride=1, padding=1, pool_win=(2,2), block=1, down_size=True,
                                     use_bn = True, use_instn = False, use_ibn = False)
        self.down_block5 = convLayer_ibn(input_channel=512, output_channel=512, filter_size=3, 
                                     stride=1, padding=1, pool_win=(2,2), block=1, down_size=True,
                                     use_bn = True, use_instn = False, use_ibn = False)
        self.down_block6 = convLayer_ibn(input_channel=512, output_channel=512, filter_size=3, 
                                     stride=1, padding=1, pool_win=(2,2), block=1, down_size=True,
                                     use_bn = True, use_instn = False, use_ibn = False)
        
        self.avgpool = torch.nn.AdaptiveAvgPool2d(1)
        
        self.flatten = torch.nn.Flatten()
        
        self.act_relu = torch.nn.ReLU()
        
        self.fc = torch.nn.Linear(32768,512)
        
        self.out = torch.nn.Linear(512,1)


    def forward(self, x):
        x1 = self.down_block1(x)
        x2 = self.down_block2(x1)
        # print(x2.shape)
        x3 = self.down_block3(x2)
        x4 = self.down_block4(x3)
        x5 = self.down_block5(x4)
        x6 = self.down_block6(x5)
        x6_1 = self.flatten(x6)
        x7 = self.fc(x6_1)
        x8 = self.act_relu(x7)
        x = self.out(x8)
        
        return x
