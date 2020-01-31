#cat reader

def cat_reader(path,first_line=False): #
    
    file=open(path,'r')
    if file.mode=='r':
        out=file.read()
        lines=out.split('\n')    
        nested_data=[line.split() for line in lines]
     #   print(nested_data)
        if first_line==True:
            output=[[float(i[j]) for i in nested_data] for j in range(len(nested_data[1]))]
          #  print(output)  
        else:
            nested_data.pop(0)
            output=[[float(i[j]) for i in nested_data] for j in range(len(nested_data[1]))]
        
        return output
        

print(cat_reader("z_puv.cat"))