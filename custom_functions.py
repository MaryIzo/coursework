import pandas as pd
import numpy as np
import scipy
import matplotlib
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
import math

from sklearn import metrics
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import LabelEncoder
from sklearn.cluster import AgglomerativeClustering
import datetime


from sklearn.cluster import DBSCAN
from sklearn.cluster import OPTICS
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import MeanShift
from sklearn.cluster import AgglomerativeClustering
import hdbscan


def matrix_of_distance(X, Y, E, label=2):
    """
    Возвращает матрицу рассстояний, посчитанную с использованием трёх разных метрик.
    label = 1 евклидово расстояние r_ij.
    label = 2 домножается на корень из энергий sqrt(e_i*e_j)*r_ij.
    label = 3 домножается на дробь (e_i*e_j)/(e_i+e_j)*r_ij.
    
    Входные значения:
    X - координата x, 
    Y - координата y,
    E - энергия, 
    label - маркер метрики расстояния.
    
    Возвращаемые значения:
    датафрейм матрицы взаимных расстояний.
    """    
    # out = попарное евклидовое расстояние
    Z = np.array([[complex(X[i], Y[i]) for i in range(len(X))]])
    out = abs(Z.T-Z)
    
    if label==1:
        # евклидовое расстояние
        oute = 1
    elif label==2:
        # sqrt(e_i*e_j)
        oute = np.sqrt(np.outer(E, E))
    elif label==3:
        # e_i*e_j/(e_i+e_j)
        outmult = np.outer(E, E)
        outadd = np.add.outer(list(E), list(E))
        oute = outmult/outadd
        
    # Для треугольной формы    
    #oute = np.triu(oute)
    return pd.DataFrame(oute*out)


def clusters_change(clusters, val_to_change, val):
    """
    Заменяет значение в словаре, который содержит номера кластеров частиц семейства
    
    Входные значения:
    clusters - словарь кластеров 
    val_to_change - заменяемое значение (не ключ) 
    val - то значение, на которое нужно заменить
    
    Выходные значения:
    clusters - словарь кластеров
    """
    for key in clusters.keys():
        if clusters[key]==val_to_change:
            clusters[key] = val
    return clusters

def df_update(dfv, min_j):
    """
    Обновляет датафрейм, удаляя из него частицу.
    
    Входные значения:
    dfv - значения датафрейма (из метода values).
    min_j - номер удаляемой частицы.
    
    Выходные значения:
    dftmp - новый датафрейм.
    """
    dftmp = pd.DataFrame(dfv)
    dftmp = dftmp[dftmp[5]!=min_j].copy().reset_index(drop=True)
    check = len(dftmp)
    dftmp[5] = [i for i in range(check)]  
    return dftmp


def PamirAlgorithm(OneFamily, 
                   eps = 48, 
                   num_name = ' num_of_family',
                   j_name = ' j',
                   x_name = 'X(J)',
                   y_name = 'Y(J)',
                   e_name = 'E(J)',
                   alg_ind = 'j_new'):
    """
    Алгоритм Памира для метрики корень из произведения энергий на евклидовое расстояние sqrt(e_i*e_j)*r_ij.
    1. Есть некоторое количество точек на плоскости, для которых известны взаимные расстояния.  
    2. Находим точки с минимальным взаимным расстоянием.  
    3. Точки объединяются в одну энергетически взвешенно  X_new = (e_i*x_i + e_j*x_j)/(e_i+e_j) и энергия e_new = (e_i+e_j).  
    4. Переходим к следующей точке, пока расстояния между точками меньше 48 или точки не закончатся. 
    Значение 48 получено в процессе работы А.С.Борисова, можно изменить для экспериментов.
    
    Входные значения:
    OneFamily - датафрейм кластеризуемого семейства, 
    eps - параметр остановки, 
    num_name - имя колонки с номерами семейств,
    j_name - имя колонуи с номерами частиц,
    x_name - имя колонки с координатой x,
    y_name - имя колонки с координатой y,
    e_name - имя колонк с энергией,
    alg_ind - имя колонки для перенумерации частиц (могут быть пропуски в значениях).
    
    Возвращаемое значение:
    labels - колонка кастеров, порядок соответствует порядку частиц в колонке с номерами частиц j_name.
    """
    
    # создадим словарь кластеров, так как нумерация частиц может быть некорректной
    points = OneFamily[j_name].unique()
    clusters = {}
    for k in points:
        clusters[k]=-1
    k = 0

    OneFamily[alg_ind] = [i for i in range(len(OneFamily))] # добавляем вспомогательный индекс для выбрасывания частиц
    dfv = OneFamily[[num_name, j_name, x_name, y_name, e_name, alg_ind]].values
    dftmp = pd.DataFrame(dfv)
    check = len(dftmp)
    
    if (check==0):
        print('Empty Family')
        return
    
    minv = np.inf
    #print(check)

    while (check!=1) and ((minv < eps) or (minv==np.inf)):

        X = dftmp[2].values
        Y = dftmp[3].values
        E = dftmp[4].values

        distm = matrix_of_distance(X, Y, E)
        l = distm.shape[0]
        distm = distm.values

        min_i = 0
        min_j = 0
        minv = np.inf
        for j in range(l):
            for i in range(j+1, l):
                if distm[j][i]< minv:
                    min_i = i
                    min_j = j
                    minv = distm[j][i]

        #print(minv, min_i, min_j)


        if minv < eps:
            # если меньше, объединяем точки
            # первая точка
            ind_i = dfv[min_i][1]
            xi = dfv[min_i][2]
            yi = dfv[min_i][3]
            ei = dfv[min_i][4]

            # вторая точка
            ind_j = dfv[min_j][1]
            xj = dfv[min_j][2]
            yj = dfv[min_j][3]
            ej = dfv[min_j][4]
            
            #print(ind_i, ind_j)
            
            if (clusters[ind_j]!=-1) and (clusters[ind_i]==-1):
                #print('old', ind_j)
                #print('max_val=', max(clusters.values()))
                clusters[ind_i] = clusters[ind_j]
                #print(clusters.values())
            elif (clusters[ind_i]!=-1) and (clusters[ind_j]==-1):
                #print('old', ind_i)
                #print('max_val=', max(clusters.values()))
                clusters[ind_j] = clusters[ind_i]
                #print(clusters.values())
            elif (clusters[ind_j]==-1) or (clusters[ind_i]==-1):
                clusters[ind_i] = k
                clusters[ind_j] = k
                k= k + 1
                #print(clusters.values())
            else:
                #print('two old')

                if (clusters[ind_i]!=clusters[ind_j]):

                    val_to_change= clusters[ind_i]
                    #print('val_to_change=', clusters[ind_i])

                    max_val = max(clusters.values())
                    #print('max_val=', max_val)

                    #print('ind_i ', clusters[ind_i], 'ind_j ', clusters[ind_j])

                    # убираем совпавший кластер
                    clusters = clusters_change(clusters, val_to_change, clusters[ind_j])

                    # заменяем максимальный номер кластера на тот, что убираем
                    clusters = clusters_change(clusters, max_val, val_to_change)
                    k = k - 1

                # после обновления значений создаём новый датафрейм
                dftmp = pd.DataFrame(dfv)
                dftmp = dftmp[dftmp[5]!=min_j].copy().reset_index(drop=True)
                check = len(dftmp)
                dftmp[5] = [i for i in range(check)]
                dfv = dftmp.values
                # проверка длины
                check = len(dftmp)
                # пересчитываем индексы

                #print(clusters.values())
                continue

            dfv[min_i][1] = ind_i # j
            dfv[min_i][2] = (ei*xi+ej*xj)/(ei+ej) # x
            dfv[min_i][3] = (ei*yi+ej*yj)/(ei+ej)  # y
            dfv[min_i][4] = (ei+ej) # e
            dfv[min_i][5] = -1  # j_new
            #print(ind_i, ind_j)
            #print(min_i, min_j)

            # после обновления значений создаём новый датафрейм
            dftmp = pd.DataFrame(dfv)
            dftmp = dftmp[dftmp[5]!=min_j].copy().reset_index(drop=True)
            check = len(dftmp)
            dftmp[5] = [i for i in range(check)]
            dfv = dftmp.values
            # проверка длины
            check = len(dftmp) 
            #print('check =', check)
        else:
            #print('minv = ', minv)
            #print('k = ', k)
            for key in clusters.keys():
                if clusters[key]==-1:
                    clusters[key] = k
                    k = k + 1
    
                    
    labels = clusters.values()
    return labels

def CoefPFI_analysis(df, num, cluster_column, hlabel_column):
    """
    Возвращает значения коэффициентов качества кластеризации в удобном для анализа формате.
    
    Входные значения:
    df - датафрейм. 
    num - номер семейства.
    cluster_column - название колонки с метками кластеров.
    hlabel_column - название колонки с метками высот.
    
    Возвращаемые значение:
    Печатает
    purity, splitting, integrity, efficiency, mean, compose
    
    (чистота, фрагментарность, целостность, эффективность, среднее, произведение средних).
    """
    
    P = 0
    E = 0
    DomH = [] #доминирующие высоты в кластерах
    N_c = len(df)
    
    labels = df[cluster_column]
    
    n_clusters_ = len(set(labels))
    
    for i in range(n_clusters_):
        
        value = df[df[cluster_column]==i][hlabel_column].value_counts()

        n_d = value[:1].values
        n_c = len(df[df[cluster_column]==i][hlabel_column])
        
        
        P_i = n_d[0]/n_c
        P = P + P_i
        
        x = value[:1].index[0]
        n_t = len(df[df[hlabel_column]==x])
        
        E_i = n_d[0]/n_t
        E = E + E_i
        DomH.append(x)
    
    purity = P/n_clusters_ #1
    
    splitting = len(set(DomH))/n_clusters_ #2
    
    integrity = len(set(DomH))/len(set(df[hlabel_column])) #3
 
    efficiency =  E/n_clusters_ #4
    
    print(f'P = {round(purity,3)}, S = {round(splitting,3)}, I = {round(integrity, 3)}, E = {round(efficiency, 3)}')
   
    
    result = [ num, round(purity, 3), round(splitting,3), round(integrity, 3), round(efficiency,3)] 
    mean = np.mean(result[1:])
    var = np.var(result[1:])
    compose = (purity+efficiency)*(integrity+splitting)/4
    
    result = result + [mean]
    result = result +[var]
    
    
    print(f'Mean: {round(mean, 3)}, compose: {round( compose , 3)}')
    
    return 


def CoefPFI(OneFamily, num, name = 'cluster'):
    """
    Возвращает значения коэффициентов качества кластеризации.
    
    Входные значения:
    OneFamily - датафрейм. 
    num - номер семейства.
    name - имя колонки с метками кластеров.
    
    Возвращаемые значение:
    result - список коэффициентов в следующем порядке:
    num, purity, splitting, integrity, efficiency, mean, var, compose
    
    (номер семейства, чистота, фрагментарность, целостность, эффективность,
    среднее, дисперсия, произведение средних).
    """
    
    P = 0
    E = 0
    DomH = [] #доминирующие высоты в кластерах
    N_c = len(OneFamily)
    
    labels = OneFamily[name]
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    
    for i in range(n_clusters_):
        value = OneFamily[OneFamily[name]==i]['labels'].value_counts()
        n_d = value[:1].values
        n_c = len(OneFamily[OneFamily[name]==i][name])
        P_i = n_d[0]/n_c
        P = P + P_i
        
        x = value[:1].index[0]
        n_t = len(OneFamily[OneFamily['labels']==x])
        E_i = n_d[0]/n_t
        E = E + E_i
        
        DomH.append(x)
    
    purity = P/n_clusters_ #1
    
    splitting = len(set(DomH))/n_clusters_ #2
    
    integrity = len(set(DomH))/len(set(OneFamily['labels'])) #3
    
    efficiency =  E/n_clusters_ #4
   
    
    result = [num, round(purity,3), round(splitting,3), round(integrity,3), round(efficiency, 3) ]
    mean = np.mean(result[1:])
    var = np.var(result[1:])
    compose = (purity+efficiency)*(integrity+splitting)/4
    
    result = result+[round(mean,3)]
    result = result + [round(var,3)]
    result = result + [round(compose,3)]
    return result


# для сокращения 

def noise_cluster(labels):
    """
    Убирает шумовые точки, обозначая их отдельными кластерами.
    
    Входные значения:
    labels - метки кластеров.
    
    Выходные значения:
    labels - метки кластеров безшума -1.
    """
    if -1 not in labels:
        return labels
    
    coef = len(set(labels))-1
    for i in range(len(labels)):  
        if labels[i] == -1:
            coef+=1
            labels[i] = labels[i]+coef
    return labels


def ClusteringAGLCLUST(eps, X):
    clustering = AgglomerativeClustering(n_clusters=None, affinity = 'precomputed', distance_threshold = eps, linkage = 'average').fit(X)
    labels = clustering.labels_    
    return labels

def ClusteringDBSCAN(eps, X):
    db = DBSCAN(eps=eps, min_samples = 1).fit(X)
    labels = db.labels_
    return noise_cluster(labels)

def ClusteringMeanShift(eps, X):
    clustering = MeanShift(bandwidth=eps).fit(X)
    labels = clustering.labels_
    return labels

def ClusteringOptics(eps, X):
    clustering = OPTICS(min_samples = 2).fit(X)
    labels = clustering.labels_
    return noise_cluster(labels)

def ClusteringHdbscan(eps, X):
    clustering = hdbscan.HDBSCAN(algorithm='best', alpha=1.0, approx_min_span_tree=True,
                            gen_min_span_tree=False, leaf_size=1, min_cluster_size=2, min_samples=None, p=None).fit(X)
    labels = clustering.labels_
    return noise_cluster(labels)

def ClusteringAffinity(eps, X):
    clustering = AffinityPropagation(damping = 0.9, random_state=5).fit(X)
    labels = clustering.labels_
    return labels

def psi(expected_array, actual_array, buckets=10):
    ''' Эта метрика позволяет измерить то, 
    насколько переменная изменилась в распределении между двумя выборками.
    Чем больше, тем сильнее изменение.
    '''
    def scale_range (array, min_arr, max_arr):
        array += -(np.min(array))
        array /= np.max(array) / (max_arr - min_arr)
        array += min_arr
        return array
    
    expected_array = np.array(expected_array)
    actual_array = np.array(actual_array)
    
    
    breakpoints = np.arange(0, buckets + 1) / (buckets) * 100
    breakpoints = scale_range(breakpoints, np.min(expected_array), np.max(expected_array))
    
    expected = np.histogram(expected_array, breakpoints)[0] / len(expected_array)
    expected = np.where(expected == 0, 0.00001, expected)
    actual = np.histogram(actual_array, breakpoints)[0] / len(actual_array)
    actual = np.where(actual == 0, 0.00001, actual)
    
    psi_value = np.sum((expected - actual)*np.log(expected/actual))
    return psi_value 
