import matplotlib.pyplot as plt
#import numpy as np
#import MPhys_model.py as mod

path=None

#PLOTS PARTICLE POSITION AS STATIC GRAPH
def plot(x,y,title,x_lab,y_lab,path=None):   
    plt.figure(figsize=(7,7))
    """
    plt.ylim(bottom=0)
    plt.ylim(top=max(y))
    """
    plt.xlim(left=0)
    plt.xlim(right=max(x)*1.1)
    plt.title(title)
    plt.xlabel(x_lab)
    plt.ylabel(y_lab)
    plt.plot(x,y)
    if path==None:
        plt.show()
    else:
        plt.savefig("plots\\"+path+"\\"+title+".png")
    

def single_plot(x,y):
    plt.plot(x,y)    
    return
    

def multiplot(x_list,y_list,title,x_lab,y_lab,legend_list=None):
    max_x=max([max(i) for i in x_list])
    
    plt.figure(figsize=(7,7))
    plt.xlim(left=0)
    plt.xlim(right=max_x*1.1)
    plt.title(title)
    plt.xlabel(x_lab)
    plt.ylabel(y_lab)
    for i in range(len(x_list)):
        single_plot(x_list[i], y_list[i])
    if legend_list != None:
        plt.legend(legend_list)
        
    if path==None:
        plt.show()
    else:
        plt.savefig("plots\\"+path+"\\"+title+".png")
    return