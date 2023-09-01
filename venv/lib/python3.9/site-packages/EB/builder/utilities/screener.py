__author__ = 'robswift'


def find_best_ensemble(results, options):
    """
    Return the best performing ensemble. If the user hasn't specified a FPF, the default behavior sorts ensembles by
    largest enrichment factor at the smallest FPF (1 / n, where n is the total number of decoys). If the user supplied
    a FPF, ensembles are sorted by the largest enrichment factor at the input FPF. If the user supplied a FPF = 1,
    ensembles are sorted by AUC. Ties are broken by considering enrichment factors at the smallest FPF not already
    considered.
    :param results: {ensemble_storage_object1, ensemble_storage_object2, ..., ensemble_storage_objectn}
    :param options: options object, contains user-specified arguments as attributes.
    :return: ensemble_storage_object (classification.EnsembleStorage) that contains the best performing ensemble
    """

    # We need the total number of decoys in the set to determine the number of decoys that correspond to the FPF values
    # at which enrichment factors were measured. The number of decoys are keys in the ef dictionary of each
    # ensemble object stored in the results dictionary.
    n = sorted(list(list(results.items())[0][1].ef.keys()), reverse=True)[0]

    # determine the number of decoys that correspond to the FPF used for training
    if not options.fpf:
        ndecoys = 1
    else:
        ndecoys = int(round(n * options.fpf))

    # sort the results according to the user-specified training method
    if ndecoys == n:
        # the user specified an fpf of 1, so wants the ensemble the maximizes the AUC, so sort on auc
        prop_key = 'auc'
        sorted_list = sorted(results.items(), key = lambda x: x[1].get_prop(prop_key), reverse=True)
    else:
        # the user is interested in an ensemble that maximizes an enrichment factor at some FPF
        prop_key = 'ef'
        sorted_list = sorted(results.items(), key = lambda x: x[1].get_prop(ndecoys, prop_key), reverse=True)

    # we only need to consider breaking a tie if there is more than one ensemble to consider
    if len(sorted_list) > 1:
        sorted_list = tie_break(sorted_list, results, prop_key, ndecoys)

    return sorted_list[0][0]


def tie_break(sorted_list, results, prop_key, ndecoys):
    """
    When finding the best ensemble, using "find_best_ensemble," in the event of a tie, we break it by looking at
    enrichment factors at successively larger FPF values, beginning with the smallest FPF not trained on.
    :param sorted_list:
    :param results:
    :return:
    """

    # Generate a list of the number of decoys that correspond to the FPF
    # values at which enrichment factors were determined, not including ndecoys, the training value
    ndecoys_list = sorted([x for x in list(results.values())[0].ef.keys() if x != ndecoys])

    # break the tie the ensemble that maximizes the enrichment factor at a given FPF (that corresponds to ndecoys)
    if prop_key == 'ef':
        index = 0
        while index < len(ndecoys_list) and sorted_list[0][1].get_prop(ndecoys_list[index], prop_key) == \
                sorted_list[1][1].get_prop(ndecoys_list[index], prop_key):
            sorted_list = sorted(results.items(), key=lambda r: (r[1].get_prop(ndecoys, prop_key),
                                                                 r[1].get_prop(ndecoys_list[index], 'ef')),
                                 reverse=True)
            index += 1
    # break the tie if the ensemble that maximizes the auc is desired
    elif prop_key == 'auc':
        index = 0
        while sorted_list[0][1].get_prop(prop_key) == sorted_list[1][1].get_prop(prop_key) \
                and index < len(ndecoys_list):
            sorted_list = sorted(results.items(), key = lambda x: (x[1].get_prop(prop_key),
                                                                   x[1].get_prop(ndecoys_list[index])))
            index += 1

    return sorted_list


def screener(molecules, ensemble, sort_order):
    """
    Uses the virtual screening scores for the receptors, or queries, specified in ensemble to sort the molecules in
    molecules in the direction specified by sort_order.
    :param molecules: a list of molecule objects (/classification/molecules.Molecules())
    :param ensemble: a tuple with receptors, or a query, that specifies a ensemble
    :param sort_order: 'asc' or 'dsc'. 'asc' sorts in ascending order (binding energy estimates) 'dsc' sorts in
    descending order (similarity scores, or binding probabilities)
    :return:
    """

    # screen
    modified_molecules = []
    for index in range(len(molecules)):
        modified_molecules.append(molecules[index])
        scores = []     #[(score, query)]

        for query in ensemble:
            scores.append( (molecules[index].GetProp(query), query))

        if sort_order == 'dsc':
            scores.sort(key=lambda x: float(x[0]), reverse=True)

        elif sort_order == 'asc':
            scores.sort(key=lambda x: float(x[0]))

        modified_molecules[index].SetProp('best_score', format(scores[0][0]))
        modified_molecules[index].SetProp('best_query', format(scores[0][1]))

    active = []
    decoy = []
    non_random = []
    for mol in modified_molecules:
        if float(mol.GetProp('best_score')) == 10000.00:
            if mol.GetProp('status') == 1:
                active.append(mol)
            else:
                decoy.append(mol)
        else:
            non_random.append(mol)

    if sort_order == 'dsc':
        non_random.sort(key=lambda mol: float(mol.GetProp('best_score')), reverse=True)
        #random.shuffle(rand)
    elif sort_order == 'asc':
        non_random.sort(key=lambda mol: float(mol.GetProp('best_score')))
        #random.shuffle(rand)
    # append the compounds with scores of 10,000 in the order active, decoy, active, ...
    rand = []
    decoy_length = len(decoy)
    active_length = len(active)
    if decoy_length > active_length:
        for a, d in zip(active, decoy[0:active_length]):
            rand.append(a)
            rand.append(d)
        for d in decoy[active_length:decoy_length]:
            rand.append(d)
    elif decoy_length < active_length:
        for a, d in zip(active[0:decoy_length], decoy):
            rand.append(a)
            rand.append(d)
        for a in active[decoy_length:active_length]:
            rand.append(a)
    elif decoy_length == active_length:
        for a, d in zip(active, decoy):
            rand.append(a)
            rand.append(d)

    modified_molecules = non_random + rand

    return modified_molecules