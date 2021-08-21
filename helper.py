#Import dependencies
#---------------------------------------------------
#
import pickle
import sys
import_all = True
if len(sys.argv) > 1:
    if str(sys.argv[1][:3]) == "pml":
        import_all = False
if import_all:
    #print("all imported")
    import numpy as np
    from matplotlib import pyplot as plt
    from copy import deepcopy
    from sympy.utilities.iterables import multiset_permutations
    from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list, to_tree, centroid, cut_tree
    from pandas import DataFrame as df
    from sklearn import decomposition
    from sklearn.cluster import KMeans
    import pandas as pd
    from multiprocessing import Pool
    from sklearn.decomposition import PCA,KernelPCA
    from sklearn.experimental import enable_iterative_imputer
    from sklearn.impute import IterativeImputer
    from sklearn.manifold import TSNE
    import threading
    from sklearn.cluster import AgglomerativeClustering
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)
    pd.set_option('display.width', None)
    pd.set_option('display.max_colwidth', None)


#CLASSES
#---------------------------------------------------
#
class Pdb:
    def __init__(self, code, type, res, chains, length, tpo_list = None, group=None, structure=None, atoms =None, residues=None, rms_cur=None, atps=[]):
        self.code = code
        self.type = type
        self.res = res
        self.chains = chains
        self.length =length
        self.file =None
        self.group = group
        #using store bio python structure feature as one of it's attributes
        self.structure = structure #biopython pdbparser.parser.get_structure
        self.atoms = atoms #undefined so far
        self.residues = residues #dictionary where keys are chains and values are residues
        self.rms_cur = rms_cur #rms distance after aligned between the align choice
        self.atps = atps
        self.tpo_list = tpo_list
    def __repr__(self, show=None):
        if show == "rms":
            return f"rms_cur: {self.rms_cur}"
        if show == "structure":
            return f"code: {self.code}; type: {self.type} res: {self.res} chains: {self.chains} length: {self.length} group: {self.group} atps: {self.atps}\nStructure: {self.structure}\nresidues: {self.residues} tpo_list: {self.tpo_list}"
        else:
            return f"code: {self.code}; type: {self.type} res: {self.res} chains: {self.chains} length: {self.length} group: {self.group} atps: {self.atps}"








#FUNCTIONS
#---------------------------------------------------
#
def nxt_itm(line,idx):
    """function which takes a line, and index the start of a series of spaces
    and returns the string right after the series of spaces before the next series of spaces"""
    k=idx
    while line[k] == " " and k<len(line):
        k+=1
    j = k
    while line[j] != " " and k<len(line):
        j+=1
    return line[k:j],j

def bin_search(arr, l, h, x):
    # Check base case
    if x>arr[h]:
        return False
    if h >= l:
        m = (h+ l) // 2
        if arr[m] == x:
            return True
        elif arr[m] > x:
            return bin_search(arr, l, m- 1, x)
        else:
            return bin_search(arr, m + 1, h, x)
    else:
        return False

###### CLUSTERING NOTEBOOK FUNCTIONS

def compute_kmeans(X,no_clusters=3):
    """Compute K means of 3 groups by default"""
    kmeans = KMeans(n_clusters=no_clusters)
    res = kmeans.fit(X)
    #print(f"res: \n {res}")
    classes = kmeans.predict(X)
    #print(f"classes: {classes}")
    return classes

def convert_classes_to_clusters(classes):
    """takes in array of classes where each index indicates the class
    outputs a cluster list which is usuable by statistics"""
    labels =list()
    out = list()
    for i,cluster in enumerate(classes):
        if cluster not in labels:
            labels.append(cluster)
            out.append([i])
        else:
            idx = labels.index(cluster)
            out[idx].append(i)
    return labels, out


def real_cluster_map(X,annotated_dict_list, chains_list):
    """plots the real cluster map and where everything should be:
    X reduced 2d matrix"""

    plt.figure(figsize=(20,20))
    colors  = ['tab:blue', 'tab:orange', 'tab:green']
    axes =X.T
    axis1 =axes[0]
    axis2 =axes[1]
    print(f"shape axis1: {np.shape(axis1)} axis2: {np.shape(axis2)}")
    for i,group in enumerate(annotated_dict_list):
        group_indices = np.array(annotated_dict_list[group])
        print(f"group: {group}")
        #group_indices = np.array([chains_list.index(item) for item in annotated_dict_list[group]])
        plt.scatter(axis1[group_indices],axis2[group_indices], s = 20, c= colors[i],label=group) 
    for i,label in enumerate(chains_list):
        plt.annotate(label, (axis1[i], axis2[i]))
    plt.legend()
    plt.title("Real Cluster Map")


def plot_cluster_map(X,misclassified,most_matched,chains_list, classes=None,cluster_labels=None,cluster_list=None, title =""):
    """Plots and colors the kmeans clusters:
    X is the coordinate matrix in this case reduced from PCA. Users either choose to fill both classes and cluster_labels, or 
    cluster list if they have it.
    Misclassified obtained by misclassification, a list of 4-tuples containing the information of misclassified conformations
    classes a list of conformations classed from 0 to n-1 classes where each classes[i] indicates which one of the n-1 classes
    the cluster is in. Cluster labels is a set of groups that exist in the cluster. Note that cluster[i] has the index of the
    original order of the clusters corresponds to most_matched[0][i] which is a list of already ordered but mapped to their 
    group names"""
    if len(most_matched[0]) > 4:
        print("Warning (to many labels): this plot kmeans function only colors for 4 clusters max")
        return
    if cluster_labels is not None and classes is not None:
        #map each cluster numerical label to the most_matched actual label:
        dictionary = {cluster_labels[i]: most_matched[0][i] for i in range(len(most_matched[0]))}
        func = lambda i: dictionary[i]
        classes_labelled = np.array(list(map(func, classes))) #list of classes/ groups the i-th conformation was classified as.
    elif cluster_list is not None:
        classes_labelled = list()
        for i in range(len(X)):
            for j,cluster in enumerate(cluster_list):
                if bin_search(cluster, 0, len(cluster)-1,i):
                    classes_labelled.append(most_matched[0][j])
                    break
        classes_labelled = np.array(classes_labelled)
        #print(f"classes_labelled: {classes_labelled}")
    else: 
        print(f'error no cluster list found')
        return
        

    ######### Plot clusters
    axes =X.T
    axis1 =axes[0]#.tolist()
    axis2 =axes[1]#.tolist()
    
    #18 colors can make 4x4 coloring max
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan','b', 'g', 'r', 'c', 'm', 'y', 'k']
    ### This can plot up
    misclassified_indices_dict_keys = list()
    for label in most_matched[0]:
        for other_label in most_matched[0]:
            if label != other_label:    
                misclassified_indices_dict_keys.append((label,other_label))
    misclassified_indices = {key: list()for key in misclassified_indices_dict_keys}
    all_misclass_conformations =list()
    for item in misclassified:
        misclassified_indices[(item[2],item[3])].append(item[0])
        all_misclass_conformations.append(item[1])
    
    plt.figure(figsize=(20,40))
    
    #plot misclassified clusters
    ax1= plt.subplot(2, 1, 1)

    for k,group in enumerate(misclassified_indices):
        ax1.scatter(axis1[misclassified_indices[group]],axis2[misclassified_indices[group]], s = 20, c= colors[k+len(most_matched[0])],label=group) 

    for i,label in enumerate(chains_list):
        if label in all_misclass_conformations:
            ax1.annotate(label, (axis1[i], axis2[i]))

    ax1.legend(title="(Predicted, Real)")
    ax1.title.set_text(f"Clustering Misclassification: {title}")

    # plot the Clusters
    ax2= plt.subplot(2, 1, 2)
    for i,item in enumerate(most_matched[0]):
        # print(f"item: {item}")
        # print(f"np.where(classes_labelled == item)[0]: {np.where(classes_labelled == item)[0]}")
        ax2.scatter(X[np.where(classes_labelled == item)[0],0], X[np.where(classes_labelled == item)[0],1], s = 20, c = colors[i],label=item)
    for i,label in enumerate(chains_list):
        if label not in all_misclass_conformations:
            #print(f"label: {label}")
            ax2.annotate(label, (axis1[i], axis2[i]))
    ax2.title.set_text(f"Clustering Results: {title} ")
    ax2.legend(title="Cluster Names")

    plt.show()


def binary_list_elimination(values):
    """Takes in a list of lists and outputs a sublist where all the lists have in common
    n^2logn comparison speed"""
    #print(f"values length: {len(values)}")
    if len(values)== 1:
        return values[0]
    if len(values) == 2:
        out = []
        #print(f"values[1]: {values[1]} values[0]: {values[0]}")
        for i in values[0]:
            if bin_search(values[1],0,len(values[1])-1,i):
                out.append(i)
        return out
    else:
        m = len(values)//2
        l1 = binary_list_elimination(values[:m])
        l2 = binary_list_elimination(values[m:])
        out = l1+ l2
        return binary_list_elimination(out)
def best_fit_elimination(cluster,annotated_dict_list):
    """takes in a single cluster list and annotated_dictionary """
    best_fit_counter_dict = {group: 0 for group in annotated_dict_list}
    for group in annotated_dict_list:
        temp =binary_list_elimination([cluster,annotated_dict_list[group]])
        #print(f"group: {group} temp: {temp}")
        best_fit_counter_dict[group] = len(temp)
    #print(f"best fit counter dict: {best_fit_counter_dict}")
    return best_fit_counter_dict

def micro_avg(matrix,labels):
    """Takes in big matrix along with it's labels and uses one-vs-rest approach to create 
    a list of mini 2x2 micro average matrix and an overall pooled statistics"""
    matrices= list()
    column_sums = np.sum(matrix,axis=0)
    row_sums = np.sum(matrix,axis=1)
    total_sum = np.sum(matrix)
    overall = np.zeros(shape=(2,2))
    for i in range(len(matrix)):
        temp =  np.zeros(shape=(2,2))
        temp[0,0] = matrix[i,i]
        temp[1,0] = column_sums[i] - matrix[i,i]
        temp[0,1] = row_sums[i] - matrix[i,i]
        temp[1,1] = total_sum - row_sums[i] - column_sums[i] + matrix[i , i]
        #print(f"temp: {temp}")
        matrices.append(df(temp.astype(int), index= [f"predicted {labels[i]}", f"predicted not {labels[i]}"], columns= [f"Real {labels[i]}", f"Real not {labels[i]}"]))
        overall = np.add(overall,temp)
    pooled_precision = overall[0,0]/ np.sum(overall, axis = 1)[0]
    pooled_recall = overall[0,0]/ np.sum(overall, axis = 0)[0]
    pooled  = df(overall.astype(int), index= ["predicted", "predicted not "], columns= ["Real", "Real not"]) 
    return matrices, pooled, pooled_precision, pooled_recall


def misclassification(chains_list,cluster_list, annotated_dict_list, most_matched,samples=10,show=False):
    """prints misclassified list based on the most_matched order, samples is amount of
    classified samples in the cluster
    RETURNS: different predicted labels, a list of lists ordered similarly to the predicted labels, 
    within this l.o.l. is a tuple: (idx,name,predicted label, real label)"""
    if show:    
        print("--"*30+ "\nMISCLASSIFCATION SAMPLES\n"+"--"*30)
    ##sampling 10 mismatches
    out = list()
    #labels = list()
    for i,label in enumerate(most_matched[0]): #look through each of the clusters
        counter = 0 
        temp = deepcopy(most_matched[0]) #remove the label from the list
        temp.remove(label)
        for conformation in cluster_list[i]: #within each of the clusters find max samples # of misclassfied examples
            if counter > samples:
                break 
            if conformation not in annotated_dict_list[label]: #conformation is actually the index. 
                #and if the conformation is not in the corresponding dict[label] value then it's misclassfied
                counter +=1
                for other_label in temp:
                    if conformation in annotated_dict_list[other_label]:
                        if show:
                            print(f"conformation: {chains_list[conformation]} predicted: {label} real: {other_label}")
                        # if label not in labels:
                        #     labels.append(label)
                        out.append((conformation,chains_list[conformation],label,other_label))
                        #else:
                        #     idx = labels.index(label)
                        #     out[idx].append((conformation,chains_list[conformation],label,other_label))
    return out

    
def statistics(cluster_list, annotated_dict_list,to_print=True):
    """Takes a list of groups unassigned, a dictionary of lists annotated
    where the keys are the groups and the lists are the supposed conformations
    and prints statistics matrix of cross classified conformations:
    where rows correspond to the predictions and columns are the actual"""
    if len(cluster_list) != len(annotated_dict_list):
        print("lengths of clusters and groups must be equal")
        return
    if to_print:
        print("--"*30+ "\nSTATISTICS\n" +"--"*30)
        for i in range(len(cluster_list)):
            print(f"cluster {i}, length: {len(cluster_list[i])}")
    a = annotated_dict_list.keys()
    all_permutations=list(multiset_permutations(a))
    fit_list =[best_fit_elimination(cluster,annotated_dict_list) for cluster in cluster_list]
    #permutation indices match cluster_list indices

    ################# PERMUTE FIT

    most_matched = (all_permutations[0],0)
    for permutation in all_permutations:
        #print(f"permutation: {permutation}")
        permutation_matches = 0
        for i in range(len(permutation)):
            #print(f"fit_list[i]: {fit_list[i]}")
            cluster_dict =fit_list[i]
            
            permutation_matches+= cluster_dict[permutation[i]]
        if most_matched[1]< permutation_matches:
            most_matched = (permutation,permutation_matches)
    #print(f"most matched: {most_matched}") #most matched tup[0]: correct permutation where index is cluster index and value is associated group

    ################  CROSS CLASSIFICATION MATRIX

    matrix = np.zeros(shape=(len(cluster_list),len(annotated_dict_list)))
    #most_matched[0] will be asignment order and we assume most_matched[0][0] be the group label for the 0 cluster_list[0]
    for i in range(len(most_matched[0])):
        for j,group2 in enumerate(most_matched[0]):
            matrix[i,j] = fit_list[i][group2]
    #print(f"matrix: {matrix}")


    ################# PRECISION, RECALL, F1

    precision_dict = {key :None for key in most_matched[0]}
    recall_dict = {key :None for key in most_matched[0]}
    
    column_sum =np.sum(matrix,axis=0)
    row_sum= np.sum(matrix, axis =1)
    #print(f"row_sum {row_sum} column_sum: {column_sum}")
    for i in range(len(matrix)):
        #print(f"most_matched[0][i]: {most_matched[0][i]}")
        precision_dict[most_matched[0][i]] = matrix[i,i]/ row_sum[i]
        recall_dict[most_matched[0][i]] = matrix[i,i]/ column_sum[i]
    F1_score_dict = {key: (2*precision_dict[key]*recall_dict[key]/(precision_dict[key] + recall_dict[key])) for key in most_matched[0]}
    p_avg = sum(precision_dict.values())/len(precision_dict)
    r_avg = sum(recall_dict.values())/len(recall_dict)
    f1_avg = sum(F1_score_dict.values())/len(F1_score_dict)
    stats_mat = list()
    for i,key in enumerate(most_matched[0]):
        stats_mat.append([precision_dict[key],recall_dict[key],F1_score_dict[key]])
        if i == len(most_matched[0])-1:
            stats_mat.append([p_avg,r_avg,f1_avg])
    #print(f"stats mat:\n {stats_mat}")
    stats_df = df(stats_mat,most_matched[0] + ["Averages"], ["Precision", "Recall", "F1 score"])
    micro_averages, pooled, pooled_precision,pooled_recall= micro_avg(matrix,most_matched[0])

    out = df(matrix.astype(int), ["Predicted "+label for label in most_matched[0]] ,["Real "+label for label in most_matched[0]])
    if to_print:
        ############## MACRO AVERAGES
        print("--"*30+ "\nMACRO AVERAGES\n" +"--"*30)
        #print(f"precision: {precision_dict}\nmacro average: {p_avg}\nrecall: {recall_dict}\nmacro average: {r_avg}\nF1 score: {F1_score_dict}\nmacro average:{f1_avg}")
        print(stats_df.to_string(max_cols=3))
        print(f"out:\n {out.to_string(max_cols=len(most_matched[0]))}")


        ############## MICRO AVERAGES

        print("--"*30+ "\nMICRO AVERAGES\n"+"--"*30)
        #print(f"len mic_avg: {len(micro_averages)}")
        for m_avg in micro_averages:
            print(f"{m_avg}")

        print(f"pooled:\n {pooled}")
        
        print(f"precision: {pooled_precision} recall: {pooled_recall}")
        print(f"F1_score: {2*pooled_precision*pooled_recall/(pooled_precision+pooled_recall)}")
    
    return most_matched,f1_avg,p_avg,r_avg


def heat_map(leaves_list,chains_list,matrix):
    """Prints heat map where rows are ordered by the clustering algorithm,
    columns are still chains list ordered"""
    rows = [chains_list[i] for i in leaves_list]
    ordered_mat = [matrix[i] for i in leaves_list]
    #print(f"rows: {rows}, chains_list {chains_list}")
    heat_frame = df(ordered_mat,rows,chains_list)
    print(f"starting heat function 2")
    #f, ax = plt.subplots(figsize=(11, 9))
    #cmap = sns.diverging_palette(230, 20, as_cmap=True)
    plt.figure(figsize=(1000,1000))
    plt.xticks(range(len(chains_list)),chains_list,rotation=90)
    plt.yticks(range(len(rows)),rows)
    plt.imshow(heat_frame, cmap='hot',interpolation="nearest")

def plot_pca(matrix,chains_list,title=""):
    """Reduces matrix to 2 dimensions using PCA and plots it"""
    np.matrix(matrix)
    pca = decomposition.PCA(n_components=2)
    reduced_matrix=pca.fit_transform(matrix)
    print(f"shape reduced matrix: {np.shape(reduced_matrix)}")
    
    axes =reduced_matrix.T

    plt.figure(figsize=(20,20))
    axis1 =axes[0].tolist()
    axis2 =axes[1].tolist()
    plt.scatter(axis1,axis2)
    for i,label in enumerate(chains_list):
        plt.annotate(label, (axis1[i], axis2[i]))
    plt.title(title)
    return reduced_matrix

def plot_tsne(matrix,chains_list,title="",perplexity = 50):
    """Reduces matrix to 2 dimensions using TSNE and plots it"""
    reduced_matrix =TSNE(n_components=2,init='pca',method='exact',perplexity=perplexity).fit_transform(matrix)
    axes =reduced_matrix.T

    plt.figure(figsize=(20,20))
    axis1 =axes[0].tolist()
    axis2 =axes[1].tolist()
    plt.scatter(axis1,axis2)
    for i,label in enumerate(chains_list):
        plt.annotate(label, (axis1[i], axis2[i]))
    plt.title(title)
    return reduced_matrix

def plot_kpca(matrix,chains_list, kernel="linear",title=""):
    """Reduces matrix to 2 dimensions using PCA and plots it
    kernel choices: {‘linear’, ‘poly’, ‘rbf’, ‘sigmoid’, ‘cosine’, ‘precomputed’}
    """
    np.matrix(matrix)
    pca = decomposition.KernelPCA(n_components=2,kernel=kernel)
    reduced_matrix=pca.fit_transform(matrix)
    print(f"shape reduced matrix: {np.shape(reduced_matrix)}")
    
    axes =reduced_matrix.T

    plt.figure(figsize=(20,20))
    axis1 =axes[0].tolist()
    axis2 =axes[1].tolist()
    plt.scatter(axis1,axis2)
    for i,label in enumerate(chains_list):
        plt.annotate(label, (axis1[i], axis2[i]))
    plt.title(title)
    return reduced_matrix




def get_list_choices(which_tree):
    """prompts user to choose trees based on chart"""
    choice_list = list()
    print(f"tree: {which_tree}")
    if which_tree == "Middle Tree":
        print("Make sure the tree is just right of the previous tree or else the statistics will be inaccurate")
    while True:
        choice = (input("left or right")).lower()
        i=0
        while choice not in ["left","right", "l","r"]:
            choice = (input("left or right")).lower()
        choice_list.append(choice)
        i+=1
        finished = input(f"done? yes?").lower()
        if finished == "y" or finished == "yes":
            break
    return choice_list

def get_tree(tree, choice_list):
    """Tree from dendrogram and list of left and right choices"""
    temp = tree
    for choice in choice_list:
        try:
            if choice in ["left", "l"]:
                temp = temp.get_left()
            else:
                temp = temp.get_right()
        except:
            print("this tree doesn't exist")
            return None
    return temp
def plot_dendrogram(title, out):
    plt.figure(figsize=(50,20))
    dn = dendrogram(out)
    plt.show()

def misclassified_vs_missing(misclassified):
    with open("num_missing_list_seg.var", "rb") as num_missing_list_seg_var:
        num_missing_list = pickle.load(num_missing_list_seg_var)
    misclassified_indices = [item[0] for item in misclassified]
    max_num_missing_list =max(num_missing_list)+1
    totals_missing = [num_missing_list.count(i) for i in range(max_num_missing_list)] #totals[i] contains # conformations with i missing residues 
    misclassified_missing = [0 for i in range(max_num_missing_list)]
    for i in misclassified_indices:
        misclassified_missing[num_missing_list[i]] +=1
    y = list()
    x = list()
    for i in range(max_num_missing_list):
        if totals_missing[i] >0: 
            y.append(misclassified_missing[i]/totals_missing[i])
            x.append(i)

    plt.figure(figsize=(20,20))
    plt.plot(x,y)
    plt.xlabel("Number of missing residues")
    plt.ylabel("Percentage of conformations misclassified")
    plt.title("Number of missing residues vs percentage misclassified")
    print(f"(#missing residues, percentage) {[(x[i],y[i]) for i in range(len(x))]}")
    print(f"miclassified_missing (list where value at i is amount of misclassified and i is amount of missing residues:\n {misclassified_missing}")
    print(f"totals_missing (list where value at i is number of conformations and i is amount of missing residues:\n {totals_missing}")



def pause_and_print(text):
    print(text)

def info(matrix,title,chains_list,annotated_dict_list,complete=True, kernel="linear",hierarchy_method = None, no_clusters=3, dist_metric = "euclidean", tsne= False, auto_cut_tree = True):
    """kernel choices: {‘linear’, ‘poly’, ‘rbf’, ‘sigmoid’, ‘cosine’}
    if hierarchy_method is None then uses kmeans plot of PCA otherwise use hierarchiacal method
    title: describes the matrix
    EG: info(matrix,"linear kernel PCA rep ward",chains_list,annotated_dict_list,kernel="linear",hierarchy_method = "ward", no_clusters=3)"""
    if complete is False:
        matrix = np.matrix(matrix)
        with open("removed_indices.var", "rb") as removed_indices_var:
            removed_indices = pickle.load(removed_indices_var)
        remaining_indices = [i for i in range(len(matrix)) if i not in removed_indices]
        matrix = matrix[np.ix_(np.array(remaining_indices,dtype=np.intp),np.array(remaining_indices,dtype=np.intp))]
        print(f"new matrix shape: {np.shape(matrix)}")
    if tsne:
        reduced_mat = plot_tsne(matrix,chains_list,title= f"T-SNE of {title} Matrix")
    else:
        reduced_mat= plot_kpca(matrix,chains_list, kernel=kernel,title=f" {kernel} Kernel of {title} Matrix")

    real_cluster_map(reduced_mat,annotated_dict_list, chains_list)
    if hierarchy_method is None:
        classes = compute_kmeans(reduced_mat,no_clusters= no_clusters)
        cluster_labels, cluster_list= convert_classes_to_clusters(classes)
        #print(f"cluster_list: {cluster_list}, cluster_labels: {cluster_labels}")
        most_matched = statistics(cluster_list,annotated_dict_list)[0]
        misclassified = misclassification(chains_list, cluster_list, annotated_dict_list, most_matched,samples = len(chains_list))
        plot_cluster_map(reduced_mat,misclassified,most_matched,chains_list,classes = classes, cluster_labels=cluster_labels, title= f"Kmeans of {title}")

    else:
        out = linkage(matrix, method = hierarchy_method, metric = dist_metric)
        thread = threading.Thread(target=plot_dendrogram,args= (title,out))
        thread.start()
        thread.join()
        print(f"Dendrogram plotted")
        
        if auto_cut_tree:
            ward_cluster = AgglomerativeClustering(linkage= hierarchy_method,n_clusters=no_clusters).fit(matrix)
            classes =ward_cluster.labels_
            cluster_labels, cluster_list = convert_classes_to_clusters(classes)
            most_matched = statistics(cluster_list,annotated_dict_list)[0]
            misclassified = misclassification(chains_list, cluster_list, annotated_dict_list, most_matched,samples = len(chains_list))
            plot_cluster_map(reduced_mat,misclassified,most_matched,chains_list,classes= classes, cluster_labels=cluster_labels, title = hierarchy_method)
        else:
            l = leaves_list(out)
            tree = to_tree(out)
            idx1, idx2 = -1 , -2
            while idx1 > idx2:
                ##############################
                #Split the trees into usable groups
                choice1 = get_list_choices("Left Most Tree")
                choice2 = get_list_choices("Middle Tree")
                tree1 = get_tree(tree, choice1)
                tree2 = get_tree(tree, choice2)
                if tree1 is not None and tree2 is not None:
                    idx1= tree1.count
                    idx2= idx1+tree2.count
            print(f"tree1: {choice1}, tree2: {choice2}")
            cluster_list = [np.sort(l[:idx1]),np.sort(l[idx1:idx2]),np.sort(l[idx2:])]
            n = tree.count
            most_matched = statistics(cluster_list,annotated_dict_list)[0]
            misclassified = misclassification(chains_list, cluster_list, annotated_dict_list, most_matched,samples =531)
            plot_cluster_map(reduced_mat,misclassified,most_matched,chains_list,cluster_list=cluster_list, title= hierarchy_method)
        misclassified_vs_missing(misclassified)


########## COORDINATE BASED MATRIX functions to help pml_script_coord.py



def reduce_seg(coords,seg_idx_list,chains_list,threshold=7):
    """Takes in the coordinates and picks out important segments
    seg_idx_list list of tuples of start and end 0 indexed
    list of chains, and threshold if a segment has more than 5 nones then the whole chain is thrown away"""
    seg_list = list() #list of conformation list of segment matrices which contains actual residue coordinates 
    #list(list(matrix()))
    for matrix in coords:
        temp = list() #list of segment matrix new matrix
        for seg in seg_idx_list:
            temp.append(matrix[seg[0]:seg[1]+1])
        seg_list.append(temp)
    removed_indices = list() #list of removed indices
    num_missing_list_seg = list() #list in the same order, each list[i] is the number of missing residues for the ith conformation
    for i,conformation in enumerate(seg_list):
        num_missing = 0
        add_conformation = True
        #print(f"i: {i}")
        for seg_mat in conformation:
            temp= sum(1 for k in seg_mat if k is None)
            #print(f"i: {i} temp: {temp}")
            num_missing += temp
            if temp >= threshold:
                if len(removed_indices) == 0 or removed_indices[-1] != i: #if the index is not already added
                    removed_indices.append(i)
                    add_conformation = False
        if add_conformation: #if not removed we keep it in the num_list_seg for statistics purposes
            num_missing_list_seg.append(num_missing)
    print(f"Threshold: {threshold} amount removed: {len(removed_indices)}")
    #print(f"segment removed_indices: {removed_indices}")
    new_seg_list = list()
    reduced_chains_list= list()
    removed_chains = list()
    if len(removed_indices)>0: #if there's any to remove
        for i,item in enumerate(seg_list):
            #print(f"i: {i} ")
            if bin_search(removed_indices,0, len(removed_indices)-1,i):
                removed_chains.append(chains_list[i])
            else:
                new_seg_list.append(item)
                reduced_chains_list.append(chains_list[i])
    else:
        reduced_chains_list = chains_list
        new_seg_list = seg_list
    with open("removed_chains.txt", "w") as removed_chains_txt:
        removed_chains_txt.write(f"Threshold: {threshold} amount removed: {len(removed_chains)}\n")
        removed_chains_txt.write(str(removed_chains))
        removed_chains_txt.close()
    with open("removed_chains.var", "wb") as removed_chains_var:
        pickle.dump([threshold]+removed_chains,removed_chains_var)
        removed_chains_var.close()
    with open("reduced_chains_list.var", "wb") as reduced_chains_list_var:
        pickle.dump(reduced_chains_list,reduced_chains_list_var)
    with open("num_missing_list_seg.var", "wb") as num_missing_list_seg_var:
        pickle.dump(num_missing_list_seg, num_missing_list_seg_var)
    with open("removed_indices.var", "wb") as removed_indices_var:
        pickle.dump(removed_indices,removed_indices_var)
    #print(f"removed indices: {removed_indices} new_seg_list: {new_seg_list} reduced_chains_list: {reduced_chains_list}")
    return new_seg_list,reduced_chains_list


def impute(matrix):
    imp = IterativeImputer()
    return imp.fit_transform(matrix)

def coordinate_impute(seg_list):
    """Impute missing values
    Parameters: seg_list: as list of conformations which are represented by a list of segment coordinate matrices, 
    (reduced by removing indices already)
    list(list(matrix()))"""
    conformation_by_seg = list() #list of segments within each segment are the coordinates 
    #for all of conformations corresponding to that segment
    for i in range(len(seg_list[0])):
        conformation_by_seg.append([conformation[i] for conformation in seg_list])
    #print(f"length conformation by seg[0]: {len(conformation_by_seg[0])}")
    conformation_by_seg_by_index = list()
    #conformation_by_seg_by_index_before_imputuation = list()
    #show = True
    imputations = 0
    for k, segment in enumerate(conformation_by_seg):
        temp = [list() for _ in range(len(segment[0]))]
        for conformation in segment:
            none_count = 0
            for i, residue in enumerate(conformation):
                if residue is None:
                    residue = np.array([np.nan, np.nan, np.nan])
                    imputations +=1
                else:
                    residue = residue[0]
                if none_count < len(conformation): #the whole array is not none
                    temp[i].append(residue)
        for i in range(len(temp)):
            #print(f"temp[i]: {temp[i]}")
            #print(f"shape temp[i]: {np.shape(temp[i])}")
            temp[i] = np.matrix(temp[i])
        #conformation_by_seg_by_index_before_imputuation.append(temp)
        #print(f"temp[-1] {temp[-1]}")
        p = Pool()
        result =  p.map(impute,temp)
        p.close()
        p.join()
        #print(f"result[-1]: {result[-1]}")
        conformation_by_seg_by_index.append(result)
    print(f"imputations done: {imputations}")
    # with open("imputed_seg_coord.var","wb") as infile:
    #     pickle.dump(conformation_by_seg_by_index,infile)
    # with open("imputed_seg_coord.txt","w") as infile2:
    #     infile2.write(str(conformation_by_seg_by_index))
    return conformation_by_seg_by_index #conformation_by_seg_by_index_before_imputuation
    #print(f"conformation by seq by index: {conformation_by_seg_by_index}")


def seg_val(conformation_by_seg_by_index,to_print=False):
    """takes seg_list and reduces it into one feature,
    coords list of matrices of coordinates
    Parameters: list(list(np.matrix()))
    A list of segment lists which contain a matrix of coordinates
    for that specific index of the sequence for each of the conformations in the same order
    for that segment
    Eg: Assume we're analysing 487 conformations of CDK2 Kinase and the segments we're analysing are
    33:44 and 150:159 respectively
    conformation_by_seg_by_index[0][0] is a matrix of coordinates of the 33rd residue 
    in each of the 487 conformations"""
    seg_list = list()
    show = True
    for i,seg in enumerate(conformation_by_seg_by_index):
        temp = [list() for _ in range(len(seg[0]))]
        for index in seg:
            for k, residue in enumerate(index):
                temp[k].append(residue)
        average = [([None]*len(seg[0][0])) for _ in range(len(seg[0]))]
        for j,conformation_seg in enumerate(temp):
            #print(f"np.sum(conformation_seg, axis=0): {np.sum(conformation_seg, axis=0)}")
            #print(f"j: {j}")
            #if j == len(temp)-1: print(f"conformation_seg: {conformation_seg} shape: {np.shape(conformation_seg)}")
            average[j] = np.sum(conformation_seg, axis=0)/(len(conformation_seg))
        if show:
            show = False
            #print(f"temp: {average}")
            #print(f"len(average): {len(average)}")
        #print(f"np.shape(average): {np.shape(average)}")
        reduced_average = [val[0] for val in TSNE(n_components=1,init="pca",perplexity = 50).fit_transform(average)]
        seg_list.append(reduced_average)
    seg_list= np.matrix(seg_list).T
    if to_print:
        print(f"seg_list: {np.shape(seg_list)}")
    with open("matrix_seg_coords.txt", "w") as matrix_seg_coords_txt:
        matrix_seg_coords_txt.write(str(seg_list))
        matrix_seg_coords_txt.close()
    with open("matrix_seg_coords.var", "wb") as matrix_seg_coords_var:
        pickle.dump(seg_list,matrix_seg_coords_var)
        matrix_seg_coords_var.close()
    return seg_list


            

def threshold_remove(threshold, segments=None):
    """Given threshold, and segments list of tuples (start and end of each segment), list(segment(start,end)),
    """
    with open("coords_res.var", "rb") as coords_res_var:
        coords = pickle.load(coords_res_var)
        coords_res_var.close()
    with open("chains_list.var", "rb") as chains_list_var:
        chains_list = pickle.load(chains_list_var)
        chains_list_var.close()
    if segments is None: #complete we add one full segment into the reduce segment analysis.
        segments = [(0,len(coords[0])-1)]
    return reduce_seg(coords,segments,chains_list,threshold=threshold)

def calc_coord_based_matrix(threshold, segments = [(34,54),(144,164)]):
    """Calculates coordinate based matrix, where each feature is the center of a pre-identified high variance region condensed by T-SNE"""
    new_seg_list, _ = threshold_remove(threshold,segments)
    #print(f"new_seg_list: {len(new_seg_list)}")
    conformation_by_seg_by_index = coordinate_impute(new_seg_list)
    return seg_val(conformation_by_seg_by_index)
    


def calc_stats_process(triple):
    """Calculate stats with different mice imputations"""
    
    new_seg_list, annotated_dict_list, no_clusters = triple
    conformation_by_seg_by_index = coordinate_impute(new_seg_list)
    #print(f"comnformation_by_seg_by_index[0][297]: {(conformation_by_seg_by_index[0][297])}")
    matrix =seg_val(conformation_by_seg_by_index)
    #thread = threading.Thread(target = pause_and_print, args= (f"matrix: {type(matrix)}"))
    #thread.start()
    #thread.join()
    ward_cluster = AgglomerativeClustering(n_clusters=no_clusters).fit(matrix)
    classes =ward_cluster.labels_
    _, cluster_list = convert_classes_to_clusters(classes)
    #print(f"shape cluster_list: {np.shape(cluster_list)}\n len(annotated_dict_list): {len(annotated_dict_list)}")
    _,f1_avg,p_avg,r_avg = statistics(cluster_list,annotated_dict_list,to_print=False)
    return [f1_avg,p_avg,r_avg]

def multirun(annotated_dict_list, segments = [(34,55),(144,165)], threshold=7,no_clusters=3, iter=10):
    """"""
    if iter <0:
        print("invalid iter")
        return
    if segments is None:
        with open("coords_res.var","rb") as coords_res_var:
            coords = pickle.load(coords_res_var)
        # for conformation in coords_res:
        #    print(f"conformation[-1]: {conformation[-1]}")
        print(f"len(coords[0])-1: {len(coords[0])-1}")
        segments = [(0,len(coords[0])-1)]
    new_seg_list, _ = threshold_remove(threshold,segments)
    vector = [(new_seg_list,annotated_dict_list,no_clusters) for _ in range(iter)]
    #print(f"vector: {vector}")
    result = list(map(calc_stats_process,vector))
    final_averages = (np.sum(np.matrix(result),axis=0)/iter)
    print(f"Iterations: {iter}, f1_avg: {final_averages[0,0]}, p_avg: {final_averages[0,1]}, r_avg: {final_averages[0,2]}")
    return result


def find_var(matrix):
    matrix = np.matrix(matrix)
    return matrix.var()


def find_variance(none_count,masked_mat):
    """Given amount of missing vectors and a masked np.nan vector, calculates the average variance of the 3 coordinates"""
    n = (len(masked_mat)-none_count)
    expectation = np.sum(masked_mat,axis=0)/n #E[X]^2
    variance = np.zeros(shape= (1,3))
    for i in range(len(masked_mat)):
        variance += np.square(masked_mat[i]-expectation)
    variance = variance/n
    return np.sum(variance)/3 # var = E[X^2] - E[X]^2, returns average var of 3 axes xyz.

def identify_highvariance_region(coords):
    """find residue to next residue of the same conformation, then align them
    together and find the variance in each residue to next residue vector
    There is a weird error with keepdims from numpy so variance is self implemented"""
    res_vect = list()
    none_count_list = list()
    for k,conformation in enumerate(coords):
        link_vectors = list()
        none_count = 0
        for i in range(len(conformation)):
            if i == len(conformation)-1:
                break
            #print(f"conformation[i][0]: {conformation[i]} conformation[i+1][0]: {(conformation[i+1])}")
            if conformation[i] is None or conformation[i+1] is None:
                link_vectors.append(np.array([0, 0, 0]))
                none_count+=1
            else:
                link_vectors.append(conformation[i+1][0]-conformation[i][0])
        res_vect.append(link_vectors)
        none_count_list.append(none_count)
    print(f"none_count_list: {none_count_list}")
    #print(f"i: {i} link_vectors[-1]: {link_vectors[-1]}")
         
    conformation_by_index_link_vectors = [list() for _ in range(len(res_vect[0]))]
    for link_vectors in res_vect:
        for i, vector in enumerate(link_vectors):
            conformation_by_index_link_vectors[i].append(vector)
    conformation_by_index_link_var = list()
    for i in range(len(conformation_by_index_link_vectors)):
        conformation_by_index_link_var.append(find_variance(none_count_list[i],conformation_by_index_link_vectors[i]))
    print(f"conformation_by_index_link_var: {(conformation_by_index_link_var)}")
    #print(f"length conformation_by_index_link_var: {len(conformation_by_index_link_var)}")
    print(f"length conformation_by_index_link_var[32:43]: {(conformation_by_index_link_var)[32:43]}")
    print(f"length conformation_by_index_link_var[149:158]: {(conformation_by_index_link_var)[149:158]}")
    return conformation_by_index_link_var


def pick_hv_regions(X,percentage):
    """Takes in the variance by index array and percentage,
    and picks all of the indices which are higher than that percentage, 
    sorted by the average variance of the subsequence"""
    mean = np.mean(X)
    std = np.sqrt(X.var())
    threshold = (1-percentage) *std + mean 
    res = list()
    for i in range(len(X)):
        if X[i] > threshold:
            res.append(i)
    print(f"res: {res}")
    subsequences = list()
    subseq= None
    for i in range(len(res)):
        if subseq is None:
            subseq = [res[i]]
        elif i+1 == len(res):
            avg = sum(X[j] for j in subseq)/len(subseq)
            subsequences.append((avg,subseq))
            break
        elif res[i] == res[i-1]+1:
            subseq.append(res[i])
        else:
            avg = sum(X[j] for j in subseq)/len(subseq)
            subsequences.append((avg,subseq))
            subseq = [res[i]]
    
    #print(f"subsequences: {subsequences}")
    subsequences.sort(key=lambda tuple: tuple[0], reverse = True)
    return subsequences


