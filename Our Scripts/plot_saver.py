import os
import time

folder="test2"

if os.path.isdir("plots\\"+folder)==False:
    os.mkdir("plots\\"+folder)

if __name__=="__main__":
    print("Generating graphs...")
    import P_L_Lya
    import EWvsZ
    import Q_Z
    print("Still generating graphs...")
    import n_ion_dot_lyc
    import t_rec
    import Q_ion_lyc
    import f_esc