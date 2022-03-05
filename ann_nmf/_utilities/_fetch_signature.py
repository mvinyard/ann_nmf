
def _fetch_signature(nmf_genes, signature):

    """
    Group and sort NMF signatures. Return the table for a single signature

    Parameters:
    -----------
    nmf_genes
        gene x signature table
        type: pandas.DataFrame

    signature
        NMF gene signature selection
        type: int or float
    """

    grouped_nmf_signatures = nmf_genes.groupby("max_id")
    return grouped_nmf_signatures.get_group(signature).sort_values(
        "max", ascending=False
    )