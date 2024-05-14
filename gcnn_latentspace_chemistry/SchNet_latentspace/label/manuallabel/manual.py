'''
Scans dataset (supports qm9) for atoms of a particular element (currently supports H,N,O) and using mol file 
and bond matrix to find the atom's connectivity and organizes into a functional group label


Inputs: 
    dataset of molecules (db file)

Outputs: 
    functional group label csv file of all atoms in selected molecules


Last update: 

    02-23-2022: no more nones! and no more unknowns! any new unkowns will now give a label of "999"
    which will automatically error out if found to let you know which molecule id has unknown 
    functional group so you can fix it in the code (or me)
    
Algorithm Summary:
    
    1) Load molecular data file (db file)
    2) Define label output file
    3) Run labelling code from DLChem package
'''

from label.manuallabel import labellerold
from label.manuallabel.utils import labellergnu, labellergnuminimal, labellerminimal, labellergnuexpand, labellergnuforacids
from schnetpack.datasets import QM9
from schnetpack import AtomsData

#load dataset
def label(qm9,dataset_filepath,label_filepath,start,end,element,gnuplot,minimal,expand,available_properties,acids):

    #run labeller on dataset
    if qm9 == True:
        qm9data = QM9(dataset_filepath,download=False,
                    remove_uncharacterized=True)
        

        if gnuplot == True:   
            if expand == True: 
                labellergnuexpand.label(qm9data,label_filepath,start,end,element) 
            if minimal == True:
                labellergnuminimal.label(qm9data,label_filepath,start,end,element)
            else:
                labellergnu.label(qm9data,label_filepath,start,end,element) 
        else:
            if minimal == True:
                labellerminimal.label(qm9data,label_filepath,start,end,element)      
            else:
                labellerold.label(qm9data,label_filepath,start,end,element) 
    
    else:

        data = AtomsData(dataset_filepath,available_properties=available_properties)
        
        if gnuplot == True: 
            if acids == True:
                labellergnuforacids.label(data,label_filepath,start,end,element)
            if expand == True: 
                labellergnuexpand.label(data,label_filepath,start,end,element) 
            if minimal == True: 
                labellergnuminimal.label(data,label_filepath,start,end,element)
            else:
                labellergnu.label(data,label_filepath,start,end,element)
        else:
            if minimal == True:
                labellerminimal.label(qm9data,label_filepath,start,end,element)  
            else:
                labellerold.label(data,label_filepath,start,end,element)
