"""Serialize / deserialize intermediate MAmotif objects with pickle."""
from __future__ import annotations

import pickle


def save_classifier(classifier, file_name):
    """Save classifier result into a file."""
    path = file_name if file_name.endswith('.classifier') else file_name + '.classifier'
    with open(path, 'wb') as handle:
        pickle.dump(classifier, handle, protocol=pickle.HIGHEST_PROTOCOL)


def read_classifier(file_path):
    """Read a classifier result from file."""
    if not file_path.endswith('.classifier'):
        print('This file is not a classifier file!')
        return None
    with open(file_path, 'rb') as handle:
        return pickle.load(handle)


def save_PeakSet(pkset, file_name):
    """Save peak set into a file."""
    path = file_name if file_name.endswith('.PeakSet') else file_name + '.PeakSet'
    with open(path, 'wb') as handle:
        pickle.dump(pkset, handle, protocol=pickle.HIGHEST_PROTOCOL)


def read_PeakSet(file_path):
    """Read PeakSet from file."""
    if not file_path.endswith('.PeakSet'):
        print('This file is not a PeakSet file!')
        return None
    with open(file_path, 'rb') as handle:
        return pickle.load(handle)


def save_GESet(geset, file_name):
    """Save gene expression set data into a file."""
    path = file_name if file_name.endswith('.GESet') else file_name + '.GESet'
    with open(path, 'wb') as handle:
        pickle.dump(geset, handle, protocol=pickle.HIGHEST_PROTOCOL)


def read_GESet(file_path):
    """Read GESet from file."""
    if not file_path.endswith('.GESet'):
        print('This file is not a GESet file!')
        return None
    with open(file_path, 'rb') as handle:
        return pickle.load(handle)


def save_GeneSet(geneset, file_name):
    """Save gene set."""
    path = file_name if file_name.endswith('.GeneSet') else file_name + '.GeneSet'
    with open(path, 'wb') as handle:
        pickle.dump(geneset, handle, protocol=pickle.HIGHEST_PROTOCOL)


def read_GeneSet(file_path):
    """Read GeneSet from file."""
    if not file_path.endswith('.GeneSet'):
        print('This file is not a GeneSet file!')
        return None
    with open(file_path, 'rb') as handle:
        return pickle.load(handle)
