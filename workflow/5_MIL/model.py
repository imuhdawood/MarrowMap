
import torch
import torch.nn as nn


class Net(nn.Module):
    def __init__(self, d):
        super(Net, self).__init__()
        self.bn1 = nn.BatchNorm1d(d,momentum=0.001)
        self.hidden1 = nn.Linear(d, d)
        self.hidden2 = nn.Linear(d//2,d//3)
        self.out = nn.Linear(d, 1)
        self.dropout = nn.Dropout(0.25)
        
    def forward(self, x):
        x = x.view(x.size(0), -1)
        # x = self.bn1(x)
        # x = torch.tanh(self.hidden1(x))
        # x = self.dropout(x)
        x = self.out(x)
        #x = torch.sigmoid(self.out(x))
        return x