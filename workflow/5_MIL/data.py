import torch
from torch.utils.data import Dataset, DataLoader

class MyDataset(Dataset):
    def __init__(self, bags):
        self.bags = bags
        
    def __getitem__(self, index):
        examples = self.bags[index]
        return torch.tensor(examples, dtype=torch.float32)

    def __len__(self):
        return len(self.bags)