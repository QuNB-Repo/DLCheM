import numpy as np
from numpy import genfromtxt
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.model_selection import train_test_split
import pandas as pd

from numpy import savetxt

from sklearn import svm
from sklearn.metrics import accuracy_score

class LDASVMwrapper:
    def __init__(self, filepath, x_cols, y_col, n_class, test_fraction,method,get_probmat=False):

        #initialize instance of variables
        self.filepath = filepath
        self.n_class = n_class
        self.x_cols = x_cols
        self.y_col = y_col
        self.test_fraction = test_fraction
        self.method = method
        self.get_probmat = get_probmat

        #Read in train and test data
        data = pd.read_csv(self.filepath)

        # Re-index the classes consecutively
        unique_labels = data['ldalabel'].unique()
        label_map = {label: i for i, label in enumerate(unique_labels)}
        data['ldalabel'] = data['ldalabel'].map(label_map)

        print(label_map)
        print(data['ldalabel'].unique())

        #process data to remove empty classes and sole members and shift class numbers
        #find classes that have less than 2 members
        class_counts = data['ldalabel'].value_counts()
        less_than_2 = class_counts[class_counts < 3].index

        #remove rows with classes that have less than 2 members
        data = data[~data['ldalabel'].isin(less_than_2)]

        # Re-index the classes consecutively
        unique_labels = data['ldalabel'].unique()
        label_map = {label: i for i, label in enumerate(unique_labels)}
        data['ldalabel'] = data['ldalabel'].map(label_map)

        print(data['ldalabel'].unique())
        print('!!!!',data['ldalabel'].value_counts()[2])
    
        # Remove rows with missing values
        data.dropna(inplace=True)

        # Define X and y
        X = data.iloc[:, 0:128]
        y = data['ldalabel']
        self.n_class = len(unique_labels)

        print(self.n_class)

        #split data, while ensuring each split gets all classes (Stratify)
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(X, y, test_size=self.test_fraction, stratify=y, random_state=42)    
        
        if self.method == 'lda':
            #Run LDA fit on training data, and save model fil in self.fit
            self.LDA = LinearDiscriminantAnalysis(store_covariance=True)
            self.fit = self.LDA.fit(self.X_train, self.y_train)

            #extract coef and intercept data of the log-posterior of the LDA, log P(k|x)
            self.coef_ = self.fit.coef_
            self.intercept_ = self.fit.intercept_
            self.covariance_ = self.fit.covariance_
            self.priors_ = self.fit.priors_

        if self.method == 'svm':
            self.svm = svm.SVC()
            self.fit = self.svm.fit(self.X_train,self.y_train)

    def FitResults(self):
        self.X_text = self.X_test.values
        self.y_test = self.y_test.values

        #calculate accuracy, extract coef and inter of fit
        self.accuracy = self.fit.score(self.X_test,self.y_test)

        #run fit on test data
        self.predicted_classes = self.fit.predict(self.X_test)

        #run testing and construct probability matrix which gives accuracy per class (on diagonal)
        if self.get_probmat == True:
            self.accuracy_perclass, self.probability_matrix = self.TestandConstructProbMat()

            if self.method == 'lda':
                return self.predicted_classes, self.accuracy, self.accuracy_perclass, self.coef_, self.intercept_, self.covariance_, self.priors_, self.probability_matrix
            if self.method == 'svm':
                return self.predicted_classes, self.accuracy, self.accuracy_perclass, self.probability_matrix
        else: 
            if self.method == 'lda':
                return self.predicted_classes, self.accuracy
            if self.method == 'svm':
                return self.predicted_classes, self.accuracy
        
    def TestResults(self,test_filepath):
        
        #load test data
        test_data = np.loadtxt(test_filepath,delimiter=',',skiprows=1)

        #parse data according to X columns and y column
        self.X_test = test_data[:, self.x_cols[0]:self.x_cols[1]]
        self.y_test = test_data[:, self.y_col]

        #run accuracy fitting score ond data
        self.accuracy = self.fit.score(self.X_test,self.y_test)

        #run prediction
        self.predicted_classes = self.fit.predict(self.X_test)
        self.predicted_classes.ke

        #run testing and construct probability matrix which gives accuracy per class (on diagonal)
        if self.get_probmat == True:
            self.accuracy_perclass, self.probability_matrix = self.TestandConstructProbMat()

            if self.method == 'lda':
                return self.predicted_classes, self.accuracy, self.accuracy_perclass, self.coef_, self.intercept_, self.covariance_, self.priors_, self.probability_matrix
            if self.method == 'svm':
                return self.predicted_classes, self.accuracy, self.accuracy_perclass, self.probability_matrix
        else: 
            if self.method == 'lda':
                return self.predicted_classes, self.accuracy
            if self.method == 'svm':
                return self.predicted_classes, self.accuracy


    def TestandConstructProbMat(self):

        #calculate accuracy per class
        self.probability_matrix = np.zeros((self.n_class,self.n_class))

        #define true classes
        true_classes = self.y_test

        #for each datapoint, calculate predicted vs true class and add 1 to the right place in probability matrix
        self.accuracy_perclass = []
        for datapoint in range(len(self.predicted_classes)):
            
            #predicted vs true class
            predicted_class = int(self.predicted_classes[datapoint])
            true_class =  int(true_classes[datapoint])

            #add one to the correct element in the probability matrix
            self.probability_matrix[predicted_class][true_class] = self.probability_matrix[predicted_class][true_class] + 1

        #normalize probability matrix, each column divide by number in class,
        #and give accuracy per class
        for _class in range(self.n_class):

            #count number of the class
            num_class = np.count_nonzero(true_classes == _class)

            if num_class != 0:

                #normalize each column of prob matrix by dividing by number in the class
                self.probability_matrix[:,_class] = np.divide(self.probability_matrix[:,_class],num_class)

                #define accuracy as the diagonal elements of probability matrix
                self.accuracy_perclass.append(self.probability_matrix[_class][_class])

            else:
                raise ZeroDivisionError('class ' + str(_class)+ ' has zero members! class')

        return  self.accuracy_perclass, self.probability_matrix

    def RootPointStudy(self, ymxb_filepath, n_features):

        ymxb = np.loadtxt(ymxb_filepath,delimiter=',')


        A = ymxb[:,0:n_features]
        b = ymxb[:,n_features]

        #need to fix b
        #invert
        cov_inv = np.linalg.inv(self.covariance_)

        #square root
        cov_inv_half = np.sqrt(cov_inv)

        #multiply by a factor 2pi^(d/2)
        twopi_cov_inv_half = cov_inv_half*(2*np.pi)**(n_features/2)
        
        #solve for root point
        X  = np.linalg.solve(A,b)

        #probability of root point
        probability = self.fit.predict_log_proba(X[:,0:n_features])


        prediction = self.fit.predict(X[:,0:n_features])
        logp = self.LDA.decision_function(X[:,0:n_features])
        
        print('prediction', prediction)
        print('incomplete logp', logp)
        print('probability of X:', probability)



