#
# Copyright (c) 2019, NVIDIA CORPORATION.  All rights reserved.
#
# NVIDIA CORPORATION and its licensors retain all intellectual property
# and proprietary rights in and to this software, related documentation
# and any modifications thereto.  Any use, reproduction, disclosure or
# distribution of this software and related documentation without an express
# license agreement from NVIDIA CORPORATION is strictly prohibited.
#

"""Contains functions to read and write to BED format."""

# Import requirements
import pandas as pd


def read_intervals(bed_file, stranded=False, skip=0):
    """Read genomic intervals from a BED file into a DataFrame.

    Args:
        bed_file: Path to BED file
        stranded: Strand information included with intervals

    Returns:
        df: Pandas DataFrame containing intervals.

    """
    if stranded:
        df = pd.read_csv(bed_file, sep='\t', header=None,
                     names=['chrom', 'start', 'end', 'strand'],
                     usecols=(0, 1, 2, 3),
                     dtype={'chrom': str, 'start': int, 'end': int,
                     'strand': str}, skiprows=skip)
    else:
        df = pd.read_csv(bed_file, sep='\t', header=None,
                     names=['chrom', 'start', 'end'],
                     usecols=(0, 1, 2),
                     dtype={'chrom': str, 'start': int, 'end': int},
                     skiprows=skip)
    return df


def read_sizes(sizes_file, as_intervals=False):
    """Read chromosome sizes into a DataFrame.

    Args:
        sizes_file(str): Path to sizes file
        as_intervals(bool): Format the DataFrame as 0-indexed intervals

    Returns:
        df: Pandas DataFrame

    """
    df = pd.read_csv(sizes_file, sep='\t', header=None, usecols=(0, 1),
                     names=['chrom', 'length'],
                     dtype={'chrom': str, 'length': int})
    if as_intervals:
        df['start'] = [0] * len(df)
        df = df[['chrom', 'start', 'length']]
        df.rename(columns={"length": "end"}, inplace=True)
    return df


def read_summits(narrowPeak_file):
    """Read summits from narrowPeak into a DataFrame.

    Args:
        narrowPeak_file: Path to narrowPeak file

    Returns:
        df: Pandas DataFrame

    """
    df = pd.read_csv(narrowPeak_file, sep='\t', header=None, usecols=(0, 1, 9),
                     names=['chrom', 'start', 'relsummit'],
                     dtype={'chrom': str, 'start': int, 'relsummit': int},
                     skiprows=1)
    df['start'] = df['start'] + df['relsummit']
    df['end'] = df['start'] + 1
    return df[['chrom', 'start', 'end']]
