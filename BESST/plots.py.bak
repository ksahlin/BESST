'''
Created on Aug 12, 2013

@author: ksahlin
'''

import os

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except ImportError, RuntimeError:
    pass


def histogram(x_, param, bins=20, x_label='x', y_label='y', title='Histogram',nr_obs = 10000):
    x = x_[:nr_obs] # If many contigs/edges we only plot 10 000 (default) for time and visuability purposes
    dest = os.path.join(param.output_directory, 'plots', title + '.png')
    plt.hist(x, bins)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)

    plt.savefig(dest)
    plt.clf()
    return()

def dot_plot(x_, y_, param, x_label='x', y_label='y', title='Dotplot', set_marker='o'):
    x = x_[:10000] # If many contigs/edges we only plot 10 000 for time and visuability purposes
    y = y_[:10000] # If many contigs/edges we only plot 10 000 for time and visuability purposes
    dest = os.path.join(param.output_directory, 'plots', title + '.png')
    plt.scatter(x, y, marker=set_marker)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    plt.grid(True)

    plt.savefig(dest)
    plt.clf()
    return


def multiple_histogram(list_of_datasets_, param, x_label='x', y_label='y', title='Stacked_histogram'):
    list_of_datasets = [list[:10000] for list in list_of_datasets_]
    dest = os.path.join(param.output_directory, 'plots', title + '.png')
    # filter out if any of the lists contains 0 elemnets
    list_of_datasets = filter(lambda x: len(x) > 0, list_of_datasets)
    for dataset in list_of_datasets:
        plt.hist(dataset, alpha=0.5)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)

    plt.savefig(dest)
    plt.clf()

    return()





#def VizualizeGraph(G, param, Information):
#    import os
#    try:
#        import matplotlib
#        matplotlib.use('Agg')
#
#        try:
#            os.mkdir(param.output_directory + '/graph_regions' + str(int(param.mean_ins_size)))
#        except OSError:
#            #directory is already created
#            pass
#        counter = 1
#        import copy
#
#        G_copy = copy.deepcopy(G)
#        RemoveIsolatedContigs(G_copy, Information)
#        CB = nx.connected_component_subgraphs(G_copy)
#        for cycle in CB:
#            nx.draw(cycle)
#            matplotlib.pyplot.savefig(param.output_directory + 'graph_regions' + str(int(param.mean_ins_size)) + '/' + str(counter) + '.png')
#            matplotlib.pyplot.clf()
#            counter += 1
#    except ImportError:
#        pass
#    return()
