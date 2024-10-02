from dscribe.descriptors import SOAP
from ase.build import molecule
from schnetpack.datasets import QM9
from numpy import savetxt



#inputs atomistic representations from ase reading of a databause, and outputs soap molecular representation of the atoms-in-molecule
def soapwrapper(qm9_true,db_filepath,start,end,savesoaps_filepath,target_element,rcut,nmax,lmax):

    #initialize SOAP
    species = ['H','C','N','O','F']

    soap = SOAP(
        species = species,
        periodic = False,
        rcut = rcut,
        nmax = nmax,
        lmax = lmax,
    )   

    #Check that you are using QM9 data base to load that one properly
    if qm9_true == True:
        qm9data = QM9(db_filepath,download=False,remove_uncharacterized=True)
    #Else, raise not yet implemented for other databases errors
    else:
        raise(NotImplementedError,"Not yet implemented for other datasets")
    

    #start an array
    soaps_outputs = []

    #go through each molecule in the dataset
    for idx in range(start,end):

        #loading bar
        if idx % 1000 == 0:
            print(idx)

        #output atomistic representations and properties
        at, props = qm9data.get_properties(idx)

        #extract number of atoms from props
        number_atoms = len(props['_atomic_numbers'])

        #check the type of element you are looking for to make comparison 
        if target_element in props['_atomic_numbers']:

            #run soap on at, atomistic system variable
            soap_molecule = soap.create(at)
            
            #go through each atom, check atomic number is the one you want, and output soap representation into a file
            for each_atom in range(number_atoms):


                #check that the element is equal to target 
                if props['_atomic_numbers'][each_atom] == target_element:
                    
                    #append soap outputs vector with soap of this atom
                    soaps_outputs.append(soap_molecule[each_atom])

    #savet the collected soaps array
    savetxt(savesoaps_filepath,soaps_outputs,delimiter=',')




