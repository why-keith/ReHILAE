"""
Saves all required plots to a unique folder
"""
import os
from datetime import datetime
import time
import multiprocessing
start=time.time()
folder="temp"

import_number=7 #number of scripts to be run

def importer_loop(i):
    #function called by multiprocessing that calls each of the required scripts in turn
    if i==0:
        import P_L_Lya
    elif i==1:
        import EWvsZ
    elif i==2:
        import Q_Z
    elif i==3:
        import n_ion_dot_lyc
    elif i==4:
        import t_rec
    elif i==5:
        import Q_ion_lyc
    elif i==6:
        import f_esc  
        
def importer(): #single thread method - depreciated
    print("Generating graphs...")
    import P_L_Lya
    import EWvsZ
    import Q_Z
   # print("Still generating graphs...")            
    import n_ion_dot_lyc
    import t_rec
    import Q_ion_lyc
    import f_esc    
    print("Time elapsed = {}s".format(time.time()-start))

if __name__=="__main__":
    
    if os.path.isdir("plots\\"+folder)==False:
        os.mkdir("plots\\"+folder)
    
    
    jobs=[]
    for i in range(import_number):
        t=multiprocessing.Process(target=importer_loop, args=(i,))#Thread(target=importer, args=(i,))
        jobs.append(t)
        t.start()
        
    for job in jobs:
        job.join()
      #  print(job)
    print("Time elapsed = {}s".format(time.time()-start))
    os.rename("plots\\"+folder,"plots\\"+str(datetime.now().strftime("%d-"+"%m-"+"%y "+"%H"+"%M"+"%S")))
