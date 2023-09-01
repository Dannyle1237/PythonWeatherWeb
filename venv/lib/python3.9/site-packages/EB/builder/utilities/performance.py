__author__ = 'robswift'

import screener

def get_sort_order(molecules):
    """
	Count up the total number of scores whose values are positve and negative.
	If a greater number are negative, then sort in ascending order (e.g. for binding energy estimates)
	Otherwise, sort in descending order (e.g. for similarity values)
	"""

    neg_count = 0
    pos_count = 0

    for index in range(len(molecules)):

        scoreList = molecules[index].GetProp('scores')
        for element in scoreList:
            if float(element) > 0:
                pos_count += 1
            elif float(element) < 0:
                neg_count += 1

    if pos_count > neg_count:
        sort_order = 'dsc'
    else:
        sort_order = 'asc'

    return sort_order


def roc_calculator(screened_molecules, status_field, active_label, decoy_label):
    """
    Calculates ROC curve
    """

    P = 0  # Total no. of actives
    N = 0  # Total no. of decoys
    tpf = [];
    tpf.append(0)  # true positive fraction list
    fpf = [];
    fpf.append(0)  # false positive fraction list
    fpindex = []  # indeces where decoys are found are labeled '1'

    # Tally the # of positives & negatives at each threshold & in total
    for index in range(len(screened_molecules)):
        if screened_molecules[index].GetProp(status_field) == active_label and index == 0:
            tpf[index] = float(1)
            P = P + 1
            fpindex.append(0)
        elif screened_molecules[index].GetProp(status_field) == active_label and index > 0:
            tpf.append(float(tpf[index - 1] + 1))
            fpf.append(float(fpf[index - 1]))
            P = P + 1
            fpindex.append(0)
        elif screened_molecules[index].GetProp(status_field) == decoy_label and index == 0:
            fpf[index] = float(1)
            N = N + 1
            fpindex.append(1)
        elif screened_molecules[index].GetProp(status_field) == decoy_label and index > 0:
            fpf.append(float(fpf[index - 1] + 1))
            tpf.append(float(tpf[index - 1]))
            N = N + 1
            fpindex.append(1)

    # calculate TPF & FPF
    for index in range(len(tpf)):
        tpf[index] = tpf[index] / P
        fpf[index] = fpf[index] / N

    return tpf, fpf, P, N


def metric_calculator(tpf, fpf, P, N, metric_list, anal=True):
    """
    Calculates VS metrics & (optionally) their 95% conficdence limits
    """

    metrics = {}

    # build dictionary to store performance metrics
    for vs_metric in [float(x) for x in metric_list if '.' in x]:
        if N >= 1 / vs_metric:
            # enrichment factor (ef) dictionary {#_of_req_decoys : [FPF,EF]}
            metrics[round(vs_metric * N)] = [vs_metric, 0]

    # non enrichment factor dictionary {'name' : ['name',value]}
    for vs_metric in [x for x in metric_list if '.' not in x]:
        metrics[vs_metric] = [vs_metric, 0, 0, 0]

    # calculate enrichment factors & auc
    metrics['AUC'][1] = 0

    for index in range(len(tpf)):
        # calculate and assign enrichment factors
        if N * fpf[index] in metrics.keys() and N * fpf[index] > 0:
            ef = tpf[index] / fpf[index]
            metrics[N * fpf[index]][1] = ef

        # calculate auc
        if fpf[index] != fpf[index - 1]:
            metrics['AUC'][1] = tpf[index] + metrics['AUC'][1]

    metrics['AUC'][1] = metrics['AUC'][1] / N

    # reshape metrics dictionary, discarding the #_of_req_decoys, i.e: {'fpf',[value,low,high]}
    reshaped = {}
    for k,v in metrics.iteritems():
            if type(v) is float or type(v) is int:
                reshaped[k] = v
            else:
                reshaped[v[0]] = v[1:]

    return reshaped


class Performance:
    """
	Raw performance. No statistical analysis
	"""

    def __init__(self):
        self.metric_list = ['0.0001', '0.001', '0.01', '0.05', 'AUC', 'Rank']

    def get_roc(self):
        return self.tpf, self.fpf

    def get_performance(self, ensemble, molecules, sort_order, ndecoys, options, anal=False):
        self.ensemble = ensemble
        self.molecules = molecules
        self.sort_order = sort_order
        self.ndecoys = ndecoys
        self.framework = options.framework
        self.status_field = options.status_field
        self.active_label = options.active_label
        self.decoy_label = options.decoy_label

        # screen molecules
        screened_molecules = screener.screener(self.molecules, self.ensemble,
                                               self.sort_order)

        # determine ROC plot & retrieve the number of actives & decoys
        self.tpf, self.fpf, P, N = roc_calculator(screened_molecules, self.status_field,
                                                  self.active_label, self.decoy_label)

        # estimate vs metrics from ROC curve
        stats = metric_calculator(self.tpf, self.fpf, P, N, self.metric_list, anal)

        # generate the ensemble ranking information
        if self.framework:
            early_unique_fw, total_fw = screener.get_frameworks(self.molecules, self.ndecoys, self.framework,
                                                                self.status_field, self.active_label, self.decoy_label)
            active_count = len(early_unique_fw)
        else:
            active_count = screener.get_active_count(self.molecules, self.ndecoys, self.status_field, self.active_label,
                                                     self.decoy_label)
        if ndecoys < N:
            stats['Rank'] = active_count
        elif ndecoys == N:
            stats['Rank'] = stats['AUC']

        return stats


