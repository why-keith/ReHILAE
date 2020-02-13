import os
from datetime import datetime

folder=(str(datetime.now().strftime("%d-"+"%m-"+"%y "+"%H"+"%M"+"%S")))

if os.path.isdir("plots\\"+folder)==False:
    os.mkdir("plots\\"+folder)

import P_L_Lya