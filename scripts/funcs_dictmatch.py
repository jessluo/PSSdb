from difflib import SequenceMatcher
from heapq import nlargest as _nlargest
import pandas as pd


def dict_match(word,iterable,n=1):
    """Use SequenceMatcher to return a list of n string(s) and index(es) of the best
    match(es) for word within the iterable object (dict, list, series).
    """

    if not n>0 :n=1
    cutoff=0.6

    result = []
    s = SequenceMatcher()
    s.set_seq2(word)
    for idx, x in enumerate(iterable):
        s.set_seq1(x)
        if s.real_quick_ratio() >= cutoff and \
           s.quick_ratio() >= cutoff and \
           s.ratio() >= cutoff:
            result.append((s.ratio(), idx,x))

    # Move the best scorers to head of list
    result = _nlargest(n, result)

    # Strip scores for the best n match(es)
    if len(result):
       return [pd.Series({"Index": index, "Words": string}) for score, index, string in result]
    else:
        return [pd.Series({"Index":-1 , "Words": ''})]