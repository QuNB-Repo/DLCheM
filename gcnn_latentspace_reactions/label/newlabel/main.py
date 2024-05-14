
from schnetpack.datasets import QM9
from schnetpack import AtomsData

from .utils import utils
import numpy as np

from numpy import savetxt


class Label():
    def __init__(self,qm9,dataset_filepath,label_filepath,start,end,elements,expand,expandinter,available_properties,colorsmarkers=[]):

        #initialize variables
        self.qm9 = qm9
        self.dataset_filepath = dataset_filepath
        self.label_filepath = label_filepath
        self.start = start
        self.end = end
        self.elements = elements
        self.expand = expand
        self.expandinter = expandinter
        self.available_properties = available_properties
        self.colorsmarkers = colorsmarkers

        #run label
        Label.RunLabel(self)
    
    def RunLabel(self):

        #get the molecule from the db file

        #if qm9 == True, use schnetpack.datasets QM9 loader, else use AtomsData loader
        if self.qm9 == True:
            data = QM9(self.dataset_filepath,download=False,remove_uncharacterized=True)
        else:
            data = AtomsData(self.dataset_filepath,available_properties=self.available_properties)
        

        
        fg_color_dictionary = {}
        fg_color_list = []
        fg_lda_dictionary = {}
        fg_lda_list = []
        fg_marker_dictionary = {}
        fg_marker_list = []

        molecule_index_list = []

        count_fg_index = 0
        #load each qm9 molecule at a time
        for each_molecule in range(self.start,self.end):
            
            #load bar
            if each_molecule % 1000 == 0:
                print(each_molecule)


            #load molecule's properties
            at, props = data.get_properties(each_molecule)

            #extract total number of atoms, atomic_number vector
            number_atoms = len(props['_atomic_numbers'])
            atomic_numbers = props['_atomic_numbers']

            #check that element chosen is in atomic_numbers
            check_element_exists = any(atomic_number in self.elements for atomic_number in atomic_numbers)

            if check_element_exists == True:
                

                #convert xyz to mol
                utils.xyz2mol(props)

                save_len_fg_color_list = len(fg_color_list)

                fg_color_dictionary, fg_color_list,fg_lda_dictionary,fg_lda_list,fg_marker_dictionary,fg_marker_list, count_fg_index = utils.labelmolfile(self.elements,self.expandinter, fg_color_dictionary, fg_color_list, fg_lda_dictionary,fg_lda_list,count_fg_index,fg_marker_dictionary,fg_marker_list)
                
                for repeat in range(save_len_fg_color_list,len(fg_color_list)):
                    molecule_index_list.append(each_molecule)
#                for each_element in self.elements:
                    #label mol file & delete temporary files
#                    utils.labelmolfile()
                    #
        print(fg_color_dictionary)
        print(fg_lda_dictionary)
        print(fg_marker_dictionary)

#        print(molecule_index_list)

        fg_label_data = np.column_stack((fg_color_list,fg_marker_list,molecule_index_list))
#        fg_label_data = np.column_stack((fg_label_data))
        savetxt('fg_label.csv',fg_label_data,delimiter=',')




#                    if atomic_numbers[each_atom] == 1:
#                        label = 0
#                        pass



        