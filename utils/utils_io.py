import os

def mkdir(dir):
    if not os.path.isdir(dir):
        os.makedirs(dir)