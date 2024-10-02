'''
These set of tools help create distributions graphs from data
'''

import numpy as np

def distributions(data,starting_point,numb_bins,increment):

#    all_dists = []

    bin = starting_point
    distribution = np.zeros((numb_bins,2))

    for each_bin in range(numb_bins):

        distribution[each_bin][0] = bin 

        for each_data in range(len(data)):

            if bin <= data[each_data][0] < bin + increment:

                distribution[each_bin][1] = distribution[each_bin][1] + 1

        bin = bin + increment

#    all_dists.append(dist[:,1])

    return distribution

def distributions_labelled(data,starting_point,numb_bins,increment):

    label = data[:,1]
    max_label = np.max(label)
    number_labels = max_label + 1

    all_distributions = []   

    for each_label in range(int(number_labels)):

        bin = starting_point
        distribution = np.zeros((numb_bins,3))
        for each_bin in range(numb_bins):

            distribution[each_bin][0] = bin 

            for each_data in range(len(data)):

                if bin <= data[each_data][0] < bin + increment and data[each_data][1] == each_label:
                    distribution[each_bin][1] = distribution[each_bin][1] + 1
                    distribution[each_bin][2] = data[each_data][3] 

            bin = bin + increment

        all_distributions.append(distribution[:,0:3])

    return all_distributions


def distributions_labelled2(data,starting_point,numb_bins,increment):

    label = data[:,1]
    max_label = np.max(label)
    number_labels = max_label + 1

    all_distributions = np.zeros((numb_bins))   

    for each_label in range(int(number_labels)):

        bin = starting_point
        distribution = np.zeros((numb_bins,3))
        for each_bin in range(numb_bins):

            distribution[each_bin][0] = bin 

            for each_data in range(len(data)):

                if bin <= data[each_data][0] < bin + increment and data[each_data][1] == each_label:
                    distribution[each_bin][1] = distribution[each_bin][1] + 1
                    distribution[each_bin][2] = data[each_data][3] 

            bin = bin + increment

        if each_label == 0:
            all_distributions = np.column_stack((all_distributions,distribution[:,0])) 
            

        all_distributions = np.column_stack((all_distributions,distribution[:,1:3]))

    return all_distributions