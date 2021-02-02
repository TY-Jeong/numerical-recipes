# Whitening: Normalization of Correlated Multivariate Data
# Following function is created with reference to the page 253-259 of text book "Neural Networks for Applied Sciences and Engineering" by Sandhya Samarasinghe


import numpy as np


def whitening(data_matrix):
    """
    Whitening - Normalization of Correlated Multivariate Data

    :param data_matrix: '(the number of data, the number of variables) size' data matrix
    :return: '(the number of data, the number of variables) size' whitened data matrix
    """
    # whitening the data

    matrix_mean = data_matrix.mean(axis=0)  # '(the number of data, the number of variables) size' data matrix
    num_data = data_matrix.shape[0]  # the number of data
    dim_data = data_matrix.shape[1]  # the number of variables

    matrix_cov = np.zeros((dim_data, dim_data))  # Covariance matrix (the covariance between the two variables)
    for i in range(num_data):
        matrix_cov = matrix_cov + np.dot((data_matrix[i, :].reshape([1, dim_data]) - matrix_mean[:].reshape([1, dim_data])).transpose(),
                                         (data_matrix[i, :].reshape([1, dim_data]) - matrix_mean[:].reshape([1, dim_data]))) / (num_data - 1.0)

    """The Covariance matrix can be transformed into a new matrix that corresponds to a set of new rescaled variables
    that have unit variance and are independent of one another. Therefore, covariance between two new variables is zero.
    This is accomplished by the well-known eigenvalue decomposition method or principal component analysis (PCA)."""

    eig_w, eig_v = np.linalg.eig(matrix_cov)
    """Each eigenvector (eig_v) defines the direction of a new axis for the rescaled variable called a principal component
    so that all new axes are perpendicular to each other. 
    Thus, this process essentially decorrelates the original input variables."""

    # To obtain the transformed variables, the original variables are multiplied by the weights (eigenvector)
    # Prior to this, weights are normalized by dividing the corresponding standart deviation (square root of eigenvalue)
    eig_v_norm = np.zeros((dim_data, dim_data))
    for i in range(data_matrix.shape[1]):
        eig_v_norm[:, i] = eig_v[:, i]/abs(eig_w[i]) ** 0.5

    # To transform the original variables to zero mean,
    # simply subtract the corresponding mean value from the values of each of the original variables.
    matrix_whitened = np.dot((data_matrix - np.dot(np.ones((num_data, 1)).reshape([num_data, 1]), matrix_mean[:].reshape([1, dim_data]))), eig_v_norm)

    return matrix_whitened
