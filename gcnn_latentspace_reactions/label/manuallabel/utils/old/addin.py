                    count_pert = 1
                    while count_pert < pert:
                        if countC != 0 or countN != 0 or countO !=0:
                            if (countC + countN) > 2:
                                for k in range(total):
                                    neighbor_index = np.array(element_indices1[k])

                                    check_if_carbon = utils.check_C(neighbor_index,name_xyz)
                                    
                                    if check_if_carbon == True:

                                        connections = utils.connection_matrix(neighbor_index,name_mol,number_atoms)
            
                                        connected_indices2 = connections[0] #since there would only be one column in the connection matrix, only one atom index is being sent in
                                    
                                        check_if_hydrogen =True
                                        while check_if_hydrogen ==True:
                                            neighbor = min(connected_indices)
                                            neighbor = int(neighbor)
                                            check_if_hydrogen = utils.check_H(neighbor,name_xyz)
                                            index = connected_indices2.index(neighbor)
                                            connected_indices2.pop(index)

                                        nn = utils.neighboring_connections(name_mol,number_atoms,neighbor2)
                                        total,countC,countH,countO,countN = utils.count_nn(nn,name_xyz)
                                        
                                        funcgroup_key = funcgrou_key + '('
                                        for i in range(countH):
                                            funcgroup_key = funcgroup_key + 'H'
                                        for i in range(countC):
                                            funcgroup_key = funcgroup_key + 'C'
                                        for i in range(countN):
                                            funcgroup_key = funcgroup_key + 'N'
                                        for i in range(countO):
                                            funcgroup_key = funcgroup_key + 'O'
                                        for i in range(countF):
                                            funcgroup_key = funcgroup_key + 'F'
                                        funcgroup_key = funcgroup_key + ')-'

                                    check_if_nitrogen = utils.check_N(neighbor_index,name_xyz)
                                    if check_if_nitrogen == True:

                                        connections = utils.connection_matrix(neighbor_index,name_mol,number_atoms)
            
                                        connected_indices2 = connections[0] #since there would only be one column in the connection matrix, only one atom index is being sent in
                                    
                                        check_if_hydrogen = True
                                        while check_if_hydrogen ==True:
                                            neighbor = min(connected_indices)
                                            neighbor = int(neighbor)
                                            check_if_hydrogen = utils.check_H(neighbor,name_xyz)
                                            index = connected_indices2.index(neighbor)
                                            connected_indices2.pop(index)

                                        nn = utils.neighboring_connections(name_mol,number_atoms,neighbor2)
                                        total,countC,countH,countO,countN = utils.count_nn(nn,name_xyz)
                                        
                                        funcgroup_key = funcgroup_key + '('
                                        for i in range(countH):
                                            funcgroup_key = funcgroup_key + 'H'
                                        for i in range(countC):
                                            funcgroup_key = funcgroup_key + 'C'
                                        for i in range(countN):
                                            funcgroup_key = funcgroup_key + 'N'
                                        for i in range(countO):
                                            funcgroup_key = funcgroup_key + 'O'
                                        for i in range(countF):
                                            funcgroup_key = funcgroup_key + 'F'
                                        funcgroup_key = funcgroup_key + ')-'
                            if (countC + countN) < 2:
                                print('hi2')
                                for k in range(len(element_indices)):
                                    print(element_indices[0])
                                    neighbor_index = np.array(element_indices[k])

                                    check_if_carbon = utils.check_C(neighbor_index,name_xyz)
                                    check_if_nitrogen = utils.check_N(neighbor_index,name_xyz)

                                    neighbor_index = [neighbor_index]

                                    if check_if_carbon == True or check_if_nitrogen == True:
                                        connections = utils.connection_matrix(neighbor_index,name_mol,number_atoms)
            
                                        connected_indices2 = connections[0] #since there would only be one column in the connection matrix, only one atom index is being sent in
            
                                        check_if_hydrogen = True
                                        while check_if_hydrogen ==True:
                                            neighbor = min(connected_indices)
                                            neighbor = int(neighbor)
                                            check_if_hydrogen = utils.check_H(neighbor,name_xyz)
                                            index = connected_indices2.index(neighbor)
                                            connected_indices2.pop(index)


                                        nn = utils.neighboring_connections(name_mol,number_atoms,neighbor)
                                        total,countC,countH,countO,countN,countF = utils.count_nn(nn,name_xyz)
                                        
                                        funcgroup_key = funcgroup_key + ''
                                        for i in range(countH):
                                            funcgroup_key = funcgroup_key + 'H'
                                        for i in range(countC):
                                            funcgroup_key = funcgroup_key + 'C'
                                        for i in range(countN):
                                            funcgroup_key = funcgroup_key + 'N'
                                        for i in range(countO):
                                            funcgroup_key = funcgroup_key + 'O'
                                        for i in range(countF):
                                            funcgroup_key = funcgroup_key + 'F'
                                        funcgroup_key = funcgroup_key + '-'


                                
                        count_pert = count_pert + 1


                            