import os

def clear(is_class=False):
    os.system("rm -f *_transfer_out*")
    os.system("rm -f *_matterpower*")
    os.system("rm -f *_conccurent_*")
    os.system("rm -f params~*")

    if is_class:
        os.system("rm -f output/*")

clear()
