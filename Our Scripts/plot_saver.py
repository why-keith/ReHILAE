import os

folder="test2"

if os.path.isdir("plots\\"+folder)==False:
    os.mkdir("plots\\"+folder)

import P_L_Lya