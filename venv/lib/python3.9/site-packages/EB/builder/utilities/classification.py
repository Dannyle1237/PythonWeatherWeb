__author__ = 'robswift'

from collections import defaultdict
from EB.builder.utilities import screener
import numpy as np
from scipy import stats


def get_sort_order(molecules):
    """
    Determine the sort order, i.e. ascending for binding energy estimates, or descending for similarity scores. First,
    count the total number of scores whose values are positive and negative. Then assign the sort order. If a greater
    number of scores are negative, an ascending sort is assigned. If a greater number of positive scores are assigned,
    a descending sort order is assigned.
    :param molecules: list [mol_object_1, mol_object_2, .... ] mol_objects are instances of common_tools.molecules
    :return: string 'asc' for ascending sort order or 'dsc' for descending sort order.
    """

    neg_count = 0
    pos_count = 0

    # count up the number of compounds with positive and negative virtual screening scores
    for index in range(len(molecules)):
        score_list = molecules[index].GetProp('scores')
        for element in score_list:
            if float(element) > 0:
                pos_count += 1
            elif float(element) < 0:
                neg_count += 1

    # compare the counts to determine the sort order
    if pos_count > neg_count:
        sort_order = 'dsc'
    else:
        sort_order = 'asc'

    return sort_order


def splitter(structured):
    """
    Separates structured data into a list of actives or a list of decoys. actives are labeled with a '1' in their status
    fields, while decoys are labeled with a '0' in their status fields.
    :param structured: either roc_structure or score_structure.
    roc_structure: list [(id, best_score, best_query, status, fpf, tpf), ..., ]
    score_structure: list [(id, best_score, best_query, status, net decoy count, net active count), ...,]
    :return: actives: list [(id, best_score, best_query, status = 1, fpf/net decoy count, tpf/net active count), ..., ]
    :return decoys: list [(id, best_score, best_query, status = 0, fpf/net decoy count, tpf/net active count), ..., ]
    """

    actives = []
    decoys = []
    for mol in structured:
        status = mol[3]
        if status == '1':
            actives.append(mol)
        elif status == '0':
            decoys.append(mol)
    return actives, decoys


def calculate_acd(acd_structure, fpf):
    """
    Determine the area under the cumulative distribution function of the active framework count, or ACD, of those
    actives ranked better than the specified False Postive Fraction.
    :param acd_structure: list [(id, best_score, best_query, status, fw smiles), ..., ]. All of the molecules in the
    list are actives.
    :param fpf: float the False Positive fraction above which the ACD is calculated
    :return: (fpf, ACD)
    """

    # initialize variables
    fw_dict = defaultdict(int) # a count of the number of molecules that share a given framework

    # loop over the molecules in acd_structure and count up the molecules consistent with each unique framework
    for mol in acd_structure:
        fw_smiles = mol[4]
        fw_dict[fw_smiles] += 1

    fw_count = [x[1] for x in fw_dict.items()] # a list of the molecule count consistent with each unique framework
    fw_count.sort(reverse=True)                # always sort from high to low s.t.: 1 >= ACD >= 0.5

    n_actives = len(acd_structure) # number of active molecules above the fpf threshold
    n_active_fw = len(fw_count)    # number of unique active frameworks above the fpf threshold
    # determine the acd
    acd = (1. / n_actives) * (1. / n_active_fw) * sum(np.cumsum(fw_count))

    return (fpf, acd)


def calculate_ef(ef_structure, fpf, diff = None, include_intervals = None):
    """
    TO DO: Add logit bounds
    Determine ROC enrichment factor (EF) and 95% confidence intervals (CI) at input false positive fraction (fpf)
    None is returned if fpf is undefined given the total no. of decoys
    :param ef_structure: list [(id, best_score, best_query, status, gamma), ..., ]
    gamma is 1 if the score of the molecule is better than the false positive threshold and 0 otherwise
    :param fpf: float
    :return (fpf, EF, EF - CI, EF + CI) or None: tuple or None
    """

    # set variables
    n = len([x for x in ef_structure if x[3] == '0'])          # no of decoys
    p = len([x for x in ef_structure if x[3] == '1'])          # no of actives


    if not fpf:
        fpf = float(1) / n

    # calculate tpf at fpf threshold tpf = E(gamma)a
    actives, decoys = splitter(ef_structure)
    tpf = 0
    for mol in actives:
        gamma = mol[4]
        tpf += gamma
    tpf = (1. / p) * tpf

    # calculate enrichment factor (ef)
    ef = tpf / fpf

    if diff:
        # active & decoy variances are independently determined in the difference calculations
        return ef, tpf
    elif include_intervals:
        # standard error
        if ef != 0:
            efvar_a, efvar_d, s2 = calculate_ef_var(tpf, fpf)
            se = (1 / fpf) * ((efvar_a / p) + s2 * (efvar_d / n) ) ** 0.5

            # 95% confidence interval
            ci = se * stats.t.ppf(1 - 0.025, n + p - 1)

            # bounds
            lower = ef - ci
            upper = ef + ci

            # dumb limit checking (ADD LOGIT FUNCTION)
            if lower < 0:
                lower = 0.0
            if upper > 1 / fpf:
                upper = 1 / fpf

            return (fpf, ef, lower, upper)
        else:
            # when the ef is 0, there's problems with taking its logarithm. This is a
            # temporary work around put in place 2/12/15.
            return(fpf, ef, 0, 0)
    else:
        decoy_threshold = int(round(n * fpf))   # number of decoys corresponding to fpf
        return (decoy_threshold, ef)


def calculate_ef_var(tpf, fpf):
    """
    determine variance due to actives (efvar_a) decoys (efvar_d) and s2, the slope of the ROC curve tangent to the
    fpf @ which the enrichment factor was calculated
    :param tpf: float tpf @ which the enrichment factor was calculated
    :param fpf: float fpf @ which the enrichment factor was calculated
    :return efvara, efvard, s2: tuple
    """

    efvara = (tpf * (1 - tpf))
    efvard = (fpf * (1 - fpf))

    ef = tpf / fpf
    if fpf == 1:
        return(0, 0, 0)
    else:
        s = ef * ( 1 + (np.log(ef)/np.log(fpf)))
        s2 = s * s
        return (efvara, efvard, s2)


def calculate_ef_diff(ef_structure1, ef_structure2, fpf):
    """
    returns the absolute value of the difference in enrichment factors and corresponding statistics.
    specifically, |ef1 - ef2|, the 95% confidence interval, and the 2-sided p-value
    :param score_structure_1: list [(id, best_score, best_query, status, net decoy count, net active count), ...,]
    :param score_structure_2: list [(id, best_score, best_query, status, net decoy count, net active count), ...,]
    :param fpf: float fpf @ which the enrichment factor was calculated
    :return EF1-EF2, EF1-EF2 - CI, EF1-EF2 + CI, p-value) : tuple
    """

    ef1, tpf1 = calculate_ef(ef_structure1, fpf, 'diff')
    ef2, tpf2 = calculate_ef(ef_structure2, fpf, 'diff')

    var1a, var1d, s2_1 = calculate_ef_var(tpf1, fpf)
    var2a, var2d, s2_2 = calculate_ef_var(tpf2, fpf)

    covara, covard = calculate_ef_covar(ef_structure1, ef_structure2)

    # standard error
    vardiffa = var2a + var1a - 2 * covara
    vardiffd = var2d * s2_2 + var1d * s2_1 + 2 * covard * (s2_2 ** 0.5) * (s2_1 ** 0.5)

    p = len([x for x in ef_structure1 if x[3] == '1'])
    n = len([x for x in ef_structure2 if x[3] == '0'])

    se = (1 / fpf) * ((vardiffa / p) + (vardiffd / n)) ** 0.5

    # 95% confidence interval
    ci = se * stats.t.ppf(1 - 0.025, n + p - 1)

    # bounds
    diff = ef1 - ef2
    lower = diff - ci
    upper = diff + ci

    # 2-sided p-value
    tt = ((diff ** 2) ** 0.5)/ se
    p_value = 2 * stats.t.sf(tt, n + p - 1)

    return (ef1, ef2, diff, (lower, upper), p_value)


def calculate_ef_covar(ef_structure1, ef_structure2):
    """
    determine the active and decoy covariance of the enrichment factors, covara, and covard, respectively.
    :param efstructure_1: list [(id, best_score, best_query, status, gamma), ...,]
    :param ef_structure2: list [(id, best_score, best_query, status, gamma), ...,]
    :return (covar_a, covar_d): tuple
    """

    # split data by activity class
    actives1, decoys1 = splitter(ef_structure1)
    actives2, decoys2 = splitter(ef_structure2)

    # covariance due to actives = E[{gamma2 - E(gamma2)a} * {gamma1 - E(gamma1)a}]a
    gamma1 = [x[4] for x in actives1]
    gamma2 = [x[4] for x in actives2]
    covara = np.cov(gamma1, gamma2)[0][1]

    # covariance due to decoys = E[{gamma2 - E(gamma2)d} * {gamma1 - E(gamma1)d}]
    gamma1 = [x[4] for x in decoys1]
    gamma2 = [x[4] for x in decoys2]
    covard = np.cov(gamma1, gamma2)[0][1]

    return covara, covard


def calculate_auc(auc_structure, sort_order, diff = None):
    """
    TO DO: Add logit bounds
    determine AUC and 95% confidence intervals (CI)
    :param roc_structure: auc_structure: list [(id, best_score, best_query, status, fpf, tpf), ..., ]
    :param sort_order: either 'asc' (binding free energy estimates) or 'dsc' (probabilities)
    :return (AUC, AUC - CI, AUC + CI): tuple
    """

    # sort by best scores
    if sort_order == 'asc':
        auc_structure.sort(key=lambda x: float(x[1]))
    elif sort_order == 'dsc':
        auc_structure.sort(key=lambda x: float(x[1]), reverse=True)

    # AUC calculation <tpf>d
    auc = 0.                         # auc average
    n = 0                           # decoy count
    for mol in auc_structure:
        status = mol[3]
        if status == '0':
            tpf = mol[5]
            n += 1
            auc += tpf
    if n == 0:
        auc = 0
    else:
        auc = (1. / n) * auc

    if diff:
        return auc
    else:
        # confidence interval calculations
        n = len([x for x in auc_structure if x[3] == '0'])          # no of decoys
        p = len([x for x in auc_structure if x[3] == '1'])          # no of actives
        aucvar_a, aucvar_d = calculate_auc_var(auc_structure)       # variance due to active and decoys
        se = (aucvar_a / p + aucvar_d / n) ** 0.5                   # standard error
        ci = se * stats.t.ppf(1 - 0.025, n + p - 1)                 # 95% confidence interval

        # bounds
        lower = auc - ci
        upper = auc + ci

        # dumb limit checking (ADD LOGIT FUNCTION)
        if lower < 0:
         lower = 0.0
        if upper > 1:
          upper = 1.0

        return (auc, lower, upper)


def calculate_auc_var(auc_structure):
    """
    determine AUC variance due to actives (aucvar_a) and decoys (aucvar_d)
    :param roc_structure: list [(id, best_score, best_query, status, fpf, tpf), ...,]
    :return (aucvar_a, aucvar_d): tuple variance due to actives and decoys, respectively
    """

    # split data by activity class
    actives, decoys = splitter(auc_structure)

    # variance due to actives = E[ (fpf - E(fpf)a) ** 2]a
    fpf = [x[4] for x in actives]
    aucvar_a = np.var(fpf, ddof=1, dtype=np.float64)

    # variance due to decoys = E[ (tpf - E(tpf)d) ** 2]d
    tpf = [x[5] for x in decoys]
    aucvar_d = np.var(tpf, ddof=1, dtype=np.float64)

    return (aucvar_a, aucvar_d)


def calculate_auc_diff(auc_structure_1, auc_structure_2, sort_order):
    """
    returns the absolute value of the difference in ROC AUC values and corresponding statistics.
    specifically, |AUC1 - AUC2|, the 95% confidence interval, and the 2-sided p-value
    :param sort_order:
    :param auc_structure_1: list [(id, best_score, best_query, status, fpf, tpf), ..., ]
    :param auc_structure_2: list [(id, best_score, best_query, status, fpf, tpf), ..., ]
    :return (AUC1-AUC2, AUC1-AUC2 - CI, AUC1-AUC2 + CI, p-value) : tuple
    """

    # determine auc and variance values for both sets
    auc1 = calculate_auc(auc_structure_1, sort_order, 'diff')
    var1a, var1d = calculate_auc_var(auc_structure_1)
    auc2 = calculate_auc(auc_structure_2, sort_order, 'diff')
    var2a, var2d = calculate_auc_var(auc_structure_2)

    # determine covariances between sets
    covara, covard = calculate_auc_covar(auc_structure_1, auc_structure_2)

    # determine standard error
    vardiffa = var2a + var1a - 2 * covara
    vardiffd = var2d + var1a - 2 * covard

    p = len([x for x in auc_structure_1 if x[3] == '1'])
    n = len([x for x in auc_structure_1 if x[3] == '0'])

    se = (vardiffa / p + vardiffd / n) ** 0.5

    # confidence interval
    ci = se * stats.t.ppf(1 - 0.025, n + p - 1)  # 95% confidence interval

    # AUC bounds
    diff = auc1 - auc2
    lower = diff - ci
    upper = diff + ci

    # 2-sided p-value prob(diff >= abs(diff)) = prob(t >= abs(tt))
    # corresponds to the null hypothesis: the two-methods perform identically
    tt = ((diff ** 2) ** 0.5) / se
    p_value = 2 * stats.t.sf(tt, n + p - 1)

    return (auc1, auc2, diff, (lower, upper), p_value)


def calculate_auc_covar(auc_structure1, auc_structure2):
    """
    determine AUC covariance due to actives (covar_a) and decoys (covar_d)
    :param auc_structure1: list [(id, best_score, best_query, status, fpf, tpf), ...,]
    :param auc_structure2: list [(id, best_score, best_query, status, fpf, tpf), ...,]
    :return (covar_a, covar_d): tuple
    """

    # split data by activity class
    actives1, decoys1 = splitter(auc_structure1)
    actives2, decoys2 = splitter(auc_structure2)

    # covariance due to actives = E[{fpf2 - E(fpf2)a} * {fpf1 - E(fpf1)a}]a
    fpf1 = [x[4] for x in actives1]
    fpf2 = [x[4] for x in actives2]
    covara = np.cov(fpf1,fpf2)[0][1]

    # covariance due to decoys = E[{tpf2 - E(tpf2)d} * {tpf1 - E(tpf1)d}]
    tpf1 = [x[5] for x in decoys1]
    tpf2 = [x[5] for x in decoys2]
    covard = np.cov(tpf1,tpf2)[0][1]    # this is only compatible with versions >= 1.5

    return covara, covard


def make_score_structure(molecules, ensemble):
    """
    puts data in the score_structure format for subsequent processing
    :param molecules: list [mol_object_1, mol_object_2, .... ] mol_objects are instances of common_tools.molecules
    :return score_structure: list [(id, best_score, best_query, status, net decoy count, net active count), ..., ]
    """
    # sort molecules by their ensemble score
    #sort_order = get_sort_order(molecules)
    sort_order = 'asc'
    sorted_molecules = screener.screener(molecules, ensemble, sort_order)

    # initiate variables
    score_structure = []
    net_active_count = 0
    net_decoy_count = 0

    for mol in sorted_molecules:
        # determine net active count & net decoy count
        status = mol.GetProp('status')
        if status == '1':
            net_active_count += 1
        elif status == '0':
            net_decoy_count += 1
        else:
            continue

        score_structure.append((mol.GetProp('id'), mol.GetProp('best_score'), mol.GetProp('best_query'), status,
                                net_decoy_count, net_active_count))

    return score_structure


def make_ef_structure(score_structure, fpf, sort_order):
    """
    Description
    :param score_structure: list [(id, best_score, best_query, status, net decoy count, net active count), ..., ]
    :param fpf: float. False positive fraction at which the enrichment factor will be calculated.
    :param sort_order: string. 'asc' (binding free energy estimate) 'dsc' (probability)
    :return: ef_structure: list [(id, best_score, best_query, status, gamma), ..., ]
    gamma is 1 if the score of the molecule is better than the false positive threshold and 0 otherwise
    """

    # set variables
    n = len([x[4] for x in score_structure if x[3] == '0']) # Total no. of decoys

    # if it wasn't defined, set the false positive fraction to the appropriate default value
    if not fpf:
        fpf = float(1) / n

    # if there are an insufficient number of decoys return None
    if n < int(round(1 / fpf)):
        return None

    # determine score threshold
    decoy_threshold = int(round(n * fpf))                                               # number of decoys corresponding
                                                                                        # to fpf
    active_list_at_threshold = [x for x in score_structure if x[4] == decoy_threshold]  # actives at ef fpf
    active_list_at_threshold.sort(key=lambda x: x[5], reverse=True)                     # sorted by active count
    score_threshold = float(active_list_at_threshold[0][1])                             # score at max active count
    ef_structure = []
    for mol in score_structure:
        score = float(mol[1])
        if sort_order == 'asc':
            if score <= score_threshold:
                gamma = 1
            elif score > score_threshold:
                gamma = 0
        elif sort_order == 'dsc':
            if score >= score_threshold:
                gamma = 1
            elif score < score_threshold:
                gamma = 0
        tup = mol[0:4] + (gamma,)
        ef_structure.append(tup)

    return ef_structure


def make_auc_structure(score_structure):
    """
    calculates ROC data points and returns roc_structure, a modified form of the score_structure format
    :param score_structure: list [(id, best_score, best_query, status, net decoy count, net active count), ..., ]
    :return auc_structure: list [(id, best_score, best_query, status, fpf, tpf), ..., ]
    """

    n = len([x for x in score_structure if x[3] == '0'])   # Total no. of decoys
    p = len([x for x in score_structure if x[3] == '1'])   # Total no. of actives

    auc_structure = []

    for mol in score_structure:
        if n == 0:
            fpf = 0
        else:
            fpf = float(mol[4]) / n
        if p == 0:
            tpf = 0
        else:
            tpf = float(mol[5]) / p
        tup = mol[0:4] + (fpf, tpf)
        auc_structure.append(tup)

    return auc_structure


def make_fpfList(options, score_structure):
    """
    aggregate default fpf values and user-defined fpf values where enrichment factor calculations will be attempted.
    :param options:
    :param score_structure:
    :return:
    """

    # set variables
    fpf = options.fpf
    n = len([x[4] for x in score_structure if x[3] == '0'])   # Total no. of decoys

    # Default fpf values. We include 1, to ensure that the ef dictionary in the results dictionary objects contains
    # the total number of decoys in the set, which is required in the screener.find_best_ensemble method
    fpf_list = [0.0001, 0.001, 0.01, 0.05, 1]

    if fpf and fpf not in fpf_list:
        fpf_list.append(float(fpf))

    if not fpf:
        fpf = float(1) / n
        fpf_list.append(fpf)

    fpf_list.sort()

    return fpf_list