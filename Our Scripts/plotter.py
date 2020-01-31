import matplotlib.pyplot as plt
#import numpy as np
#import MPhys_model.py as mod

def plot(x,y,title,x_lab,y_lab): #plots particle position as a static graph    
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
    plt.show()
    
    