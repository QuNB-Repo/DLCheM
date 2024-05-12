
def countatomsinmol(mol):

    numberatoms = 0 
    for line_index, each_line in enumerate(mol.split('\n')):
        if  line_index > 3:
            if ' H ' in each_line or ' C ' in each_line  or  ' O ' in each_line  or ' N ' in each_line  or ' F ' in each_line:
                if 'END' not in each_line:
                    numberatoms = numberatoms + 1 
    
    numberbonds = 0
    for line_index, each_line in enumerate(mol.split('\n')):
        if line_index > 3 + numberatoms:
            if 'END' not in each_line:
                numberbonds = numberbonds + 1
            if 'END' in each_line:
                break


    return numberatoms,numberbonds




#CAAAAN DELETE BELOW LATERRR!!!

def countatomsinmol(mol):

    numberatoms = 0 
    for line_index, each_line in enumerate(mol.split('\n')):
        if  line_index > 3:
            if ' H ' in each_line or ' C ' in each_line  or  ' O ' in each_line  or ' N ' in each_line  or ' F ' in each_line:
                if 'END' not in each_line:
                    numberatoms = numberatoms + 1 
    
    numberbonds = 0
    for line_index, each_line in enumerate(mol.split('\n')):
        if line_index > 3 + numberatoms:
            if 'END' not in each_line:
                numberbonds = numberbonds + 1
            if 'END' in each_line:
                break


    return numberatoms,numberbonds

def removeenv(tempmolfilepath,numberatoms,atomindex,ndepth):

    mol = open(tempmolfilepath,mode='r')
    
    mol, nbridxs = deleteatomconnections(mol,atomindex,numberatoms)
    numberatoms, numberbonds = countatomsinmol(mol)

    nbridxs2 = []
    for each_depth in range(ndepth+1):
        if each_depth == 1: 
            for each_nbr in range(len(nbridxs)):
                mol, nbrsofthisnbr = deleteatomconnections(mol,nbridxs[each_nbr],numberatoms)
                numberatoms, numberbonds = countatomsinmol(mol)

                for each_nbrnrbr in range(len(nbrsofthisnbr)):
                    nbridxs2.append(nbrsofthisnbr[each_nbrnrbr])


        if each_depth > 1:
            nbridxs = nbridxs2
            nbridxs2.clear() 
            for each_nbr in range(len(nbridxs)):

                mol, nbrsofthisnbr = deleteatomconnections(mol,nbridxs[each_nbr],numberatoms)
                numberatoms, numberbonds = countatomsinmol(mol)

                for each_nbrnrbr in range(len(nbrsofthisnbr)):
                    nbridxs2.append(nbrsofthisnbr[each_nbrnrbr])
    
    #make final adjustments
    final_mol = ''
    for line_index, each_line in enumerate(mol.split('\n')):
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
            

def deleteatomconnections(mol,atomindex,numberatoms):
    nbridxs = []
    
    #find connected atoms, delete connections, 
    mol2 = ''
    for line_index, each_line in enumerate(mol):

        if line_index > 3 + numberatoms:
            
            if ' ' +str(atomindex)+' ' in each_line[0:7]:

                if ' ' +str(atomindex)+' ' in each_line[0:4]:
                    nbridxs.append(int(each_line[4:7]))

                if ' ' +str(atomindex)+' ' in each_line[4:7]:
                    nbridxs.append(int(each_line[0:4]))
                
                each_line = each_line.replace(each_line,'')
        mol2 = mol2 + each_line
    
    #delete connected atoms
    mol3 = ''
    for each_nbr in range(len(nbridxs)):
        print(nbridxs[each_nbr])
        for line_index, each_line in enumerate(mol2):

            if  3 + numberatoms >= line_index > 3:
                if line_index == nbridxs[each_nbr] + 3:
                    each_line = each_line.replace(each_line,'')

            mol3 = mol3 + each_line


    return mol3, nbridxs
