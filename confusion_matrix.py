import numpy as np


def plot_confusion_matrix(dataset_name, cm,
                          target_names=None,
                          title='Confusion matrix',
                          cmap=None,
                          normalize=True,
                          val_auc=None):
    """
    given a sklearn confusion matrix (cm), make a nice plot

    Arguments
    ---------
    dataset_name: name of dataset
    cm:           confusion matrix from sklearn.metrics.confusion_matrix

    target_names: given classification classes such as [0, 1, 2]
                  the class names, for example: ['high', 'medium', 'low']

    title:        the text to display at the top of the matrix

    cmap:         the gradient of the values displayed from matplotlib.pyplot.cm
                  see http://matplotlib.org/examples/color/colormaps_reference.html
                  plt.get_cmap('jet') or plt.cm.Blues

    normalize:    If False, plot the raw numbers
                  If True, plot the proportions
    
    val_auc:      Value of Area Under Curve

    Usage
    -----
    plot_confusion_matrix(cm           = cm,                  # confusion matrix created by
                                                              # sklearn.metrics.confusion_matrix
                          normalize    = True,                # show proportions
                          target_names = y_labels_vals,       # list of names of the classes
                          title        = best_estimator_name) # title of graph

    Citiation
    ---------
    http://scikit-learn.org/stable/auto_examples/model_selection/plot_confusion_matrix.html

    """
    import matplotlib.pyplot as plt
    import numpy as np
    import itertools

    val_acc = np.trace(cm) / float(np.sum(cm))
    val_loss = 1 - val_acc

    if normalize:
        cm = cm.astype('float') * 100.0/ cm.sum(axis=1)[:, np.newaxis]
    
    if cmap is None:
        cmap = plt.get_cmap('Blues')

    plt.figure(figsize=(4, 2))
    plt.imshow(cm, interpolation='nearest', cmap=cmap)
    plt.title(title, pad=14)
    plt.colorbar()
    if target_names is not None:
        tick_marks = np.arange(len(target_names))
        plt.xticks(tick_marks, target_names, rotation=45)
        plt.yticks(tick_marks, target_names)

    thresh = cm.max() / 1.5 if normalize else cm.max() / 2
    for i, j in itertools.product(range(cm.shape[0]), range(cm.shape[1])):
        if normalize:
            plt.text(j, i, "{:0.2f}%".format(cm[i, j]),
                     horizontalalignment="center",
                     color="white" if cm[i, j] > thresh else "black")
        else:
            plt.text(j, i, "{:,}".format(cm[i, j]),
                     horizontalalignment="center",
                     color="white" if cm[i, j] > thresh else "black")


    plt.tight_layout()
    plt.ylabel('Predicted label')
    if normalize:
        plt.xlabel('True label\nval_acc={:0.2f}%;val_auc={:0.2f}%\nUnit: Percent of Accuracy'.format(val_acc*100.0, val_auc*100.0))
    else:
        plt.xlabel('True label\nval_acc={:0.3f};val_auc={:0.3f}\nUnit: Average accuracy in percentage'.format(val_acc, val_auc))
    # plt.savefig(dataset_name + '.png')
    plt.show()

if __name__ == '__main__':
    v_data = [
        {
            'name': 'Cirgene',
            'conf': [[ 10.11, 1.29,],[  1.69, 10.11]],
            'auc': 0.941
        },
        {
            'name': 'Colgene',
            'conf': [[ 3.44, 0.93,],[  1.36, 6.37]],
            'auc': 0.865
        },
        {
            'name': 'Ibdgene',
            'conf': [[ 1.29, 0.42,],[  1.21, 8.08]],
            'auc': 0.895
        },
        {
            'name': 'Obegene',
            'conf': [[ 14.71, 7.53,],[  1.69, 1.37]],
            'auc': 0.576
        },
        {
            'name': 'T2dgene',
            'conf': [[ 8.34, 4.76,],[  8.66, 12.64]],
            'auc': 0.673
        },
        {
            'name': 'Wt2dgene',
            'conf': [[ 3.54, 1.65,],[  1.76, 2.65]],
            'auc': 0.692
        }
    ]
    for dataset in v_data:
        plot_confusion_matrix(
                            dataset_name = dataset['name'],
                            cm           = np.array(dataset['conf']), 
                            normalize    = True,
                            target_names = ['patient', 'control'],
                            title        = "Average Confusion matrix of " + dataset['name'],
                            val_auc      = dataset['auc'])