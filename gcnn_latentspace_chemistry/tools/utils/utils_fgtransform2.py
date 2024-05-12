#This code contains the functions that perform functional group/ chemical env removal to calculate signals of those removal per functional group depth label


def count_atoms_and_bonds_mol(mol):

    numberatoms = 0 
    for line_index, each_line in enumerate(mol.split('\n')):
        if  line_index > 2: #in this case 2, because a space has been removed
            if ' H ' in each_line or ' C ' in each_line  or  ' O ' in each_line  or ' N ' in each_line  or ' F ' in each_line:
                if 'END' not in each_line:
                    numberatoms = numberatoms + 1 
    
    numberbonds = 0
    for line_index, each_line in enumerate(mol.split('\n')):
        if line_index > 2 + numberatoms: # in this case 2, ...  " " ... " " 
            if 'END' not in each_line:
                numberbonds = numberbonds + 1
            if 'END' in each_line:
                break

    return numberatoms,numberbonds

def adjust_atoms_bonds_mol(mol,numberatoms,numberbonds):
    #make final adjustments
    adjust_spaced_mol = ''
    for line_index, each_line in enumerate(mol.split('\n')):
        if line_index == 1:
            each_line = each_line.replace(each_line,each_line+'\n')
        adjust_spaced_mol = adjust_spaced_mol + each_line + '\n'
    
    final_mol = ''
    for line_index, each_line in enumerate(adjust_spaced_mol.split('\n')):
        if line_index == 3:
            if numberatoms < 10 and numberbonds < 10:
                each_line = each_line.replace(each_line[0:6],'  '+str(numberatoms)+'  '+str(numberbonds)+'')
            if numberatoms > 9 and numberbonds < 10:
                each_line = each_line.replace(each_line[0:6],' '+str(numberatoms)+'  '+str(numberbonds)+'')
            if numberatoms < 10 and numberbonds > 9:
                each_line = each_line.replace(each_line[0:6],'  '+str(numberatoms)+' '+str(numberbonds)+'')
            if numberatoms > 9 and numberbonds > 9:
                each_line = each_line.replace(each_line[0:6],' '+str(numberatoms)+' '+str(numberbonds)+'')
        final_mol = final_mol + each_line + '\n'

    return final_mol

def remove_empty_lines(input_string):
    # Split the input string into lines, filter out empty lines, and join the non-empty lines back together
    lines = input_string.split('\n')
    non_empty_lines = [line for line in lines if line.strip()]
    return '\n'.join(non_empty_lines)


def delete_atom_connections(mol,targetidx,tag):

    new_mol = ''
    new_nbridxs = []
    for line_index, each_line in enumerate(mol.split('\n')):
        
        # THE TARGET ATOM SO THAT LATER YOU CAN USE THE  TO PUT THE TARGET
        #TARGET ATOM FIRST IN THE XYZ... SO THAT IT CAN BE FOUND IN EMBEDDING DIFF
        if tag == True:
            if line_index == 3 + targetidx:
                each_line = each_line.replace(each_line,'TAG ' + each_line)

        if line_index > 3 + targetidx:
            if str(targetidx) in each_line[0:7]:

                if ' ' +str(targetidx)+' ' in each_line[0:4]:
                    new_nbridxs.append(int(each_line[4:7]))
                    each_line = each_line.replace(each_line,'')

                if ' ' +str(targetidx)+' ' in each_line[4:7]:
                    new_nbridxs.append(int(each_line[0:4]))
                    each_line = each_line.replace(each_line,'')

        new_mol = new_mol + each_line + '\n'

    return new_mol, new_nbridxs

def delete_connected_atom_and_connections(mol,targetidx):
    new_mol = ''
    for line_index, each_line in enumerate(mol.split('\n')):
        if line_index == 3 + targetidx:
            each_line = each_line.replace(each_line,'')

        new_mol = new_mol + each_line + '\n'

    new_mol, new_nbridxs = delete_atom_connections(new_mol,targetidx,tag=False)

    return new_mol, new_nbridxs


def convert_mol2xyz(mol,numberatoms):

    #First find the tag and write that first

    xyz = str(numberatoms) + '\n0.0000 \n'
    for line_index, each_line in enumerate(mol.split('\n')):
        if 'TAG' in each_line:
            xyz = xyz + each_line[35:37] + each_line[5:35]  + '\n'

    for line_index, each_line in enumerate(mol.split('\n')):
        if 3 + numberatoms >= line_index > 3 and 'TAG' not in each_line:
            xyz = xyz + each_line[31:32] + each_line[0:30] + '\n'


    return xyz


                                                    

def remove_env(tempmolfilepath,numberatoms,atomidx,ndepths):

    #offset
    ndepths = ndepths + 1

    mol = open(tempmolfilepath,mode='r')
    init_mol_string = ''
    for line_index, each_line in enumerate(mol):
        init_mol_string  = init_mol_string + each_line

    new_mol = ''
    for each_depth in range(ndepths):

        if each_depth == 0:
            #the first "neighbor" will be the target atom index
            nbridxs = [atomidx]
            for nbridxs_index, each_nbr_index in enumerate(nbridxs):

                new_mol, nbridxs2 = delete_atom_connections(init_mol_string,each_nbr_index,tag=True)


        if each_depth > 0:

            nbridxs = nbridxs2
            nbridxs2 = []
            for nbridxs_index, each_nbr_index in enumerate(nbridxs):

                new_mol, nbrnbr = delete_connected_atom_and_connections(new_mol,each_nbr_index)
                
                for nbrnbr_index, each_nbr in enumerate(nbrnbr): 
                    nbridxs2.append(each_nbr)

    new_mol = remove_empty_lines(new_mol)
    new_numberatoms,new_numberbonds = count_atoms_and_bonds_mol(new_mol)
    final_mol = adjust_atoms_bonds_mol(new_mol,new_numberatoms,new_numberbonds)
    
    final_xyz = convert_mol2xyz(final_mol,new_numberatoms)
    return final_xyz



    

