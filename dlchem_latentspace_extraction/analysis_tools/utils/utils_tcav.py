import numpy as np
from numpy import genfromtxt
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.model_selection import train_test_split

class LDAwrapper:
    def __init__(self, filepath, x_cols, y_col, n_class, test_fraction):

        #initialize instance of variables
        self.filepath = filepath
        self.n_class = n_class
        self.x_cols = x_cols
        self.y_col = y_col
        self.test_fraction = test_fraction

        #Read in train and test data
        train_data = np.loadtxt(self.filepath, delimiter=',',skiprows=1)

        #define X and Y
        X = train_data[:, self.x_cols[0]:self.x_cols[1]]
        y = train_data[:, self.y_col]

        #split data, while ensuring each split gets all classes (Stratify)
        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(X, y, test_size=self.test_fraction, stratify=y, random_state=30)    
        
        #Run LDA fit on training data, and save model fil in self.fit
        self.LDA = LinearDiscriminantAnalysis(store_covariance=True)
        self.fit = self.LDA.fit(self.X_train, self.y_train)

        #extract coef and intercept data of the log-posterior of the LDA, log P(k|x)
        self.coef_ = self.fit.coef_
        self.intercept_ = self.fit.intercept_
        self.covariance_ = self.fit.covariance_
        self.priors_ = self.fit.priors_

        return

    def FitResults(self):
        #calculate accuracy, extract coef and inter of fit
        self.accuracy = self.fit.score(self.X_test,self.y_test)

        self.predicted_classes = self.fit.predict(self.X_test)
        #run testing and construct probability matrix which gives accuracy per class (on diagonal)
#        self.accuracy_perclass, self.probability_matrix = self.TestandConstructProbMat()

        return self.accuracy, self.predicted_classes
    
    def TestResults(self,test_filepath):
        
        #load test data
        test_data = np.loadtxt(test_filepath,delimiter=',',skiprows=1)

        #parse data according to X columns and y column
        self.X_test = test_data[:, self.x_cols[0]:self.x_cols[1]]
        self.y_test = test_data[:, self.y_col]

        #run accuracy fitting score ond data
#        self.accuracy = self.fit.score(self.X_test,self.y_test)

        #run fit on test data
        self.predicted_classes = self.fit.predict(self.X_test)

        #run test and construct probability matrix of test data, which gives accuracy on diagonals
#        self.accuracy_perclass, self.probability_matrix = self.TestandConstructProbMat()

        return self.predicted_classes


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




