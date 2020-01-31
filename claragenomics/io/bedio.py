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
                     'strand': str}, skiprows=args.skip)
    else:
        df = pd.read_csv(bed_file, sep='\t', header=None,
                     names=['chrom', 'start', 'end'],
                     usecols=(0, 1, 2),
                     dtype={'chrom': str, 'start': int, 'end': int},
                     skiprows=args.skip)
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
                     dtype={0: str, 1: int})
    if as_intervals:
        df[2] = [0] * len(df)
        df.rename(columns={0: 0, 2: 1, 1: 2}, inplace=True)
        df.columns = ['chrom', 'start', 'end']
    else:
        df.columns = ['chrom', 'length']
    return df


def read_summits(narrowPeak_file):
    """Read summits from narrowPeak into a DataFrame.

    Args:
        narrowPeak_file: Path to narrowPeak file

    Returns:
        df: Pandas DataFrame

    """
    df = pd.read_csv(narrowPeak_file, sep='\t', header=None, usecols=(0, 1, 9),
                     dtype={0: str, 1: int, 9: int},
                     columns={'chrom', 'start', 'relsummit'})
    df['start'] = df['start'] + df['relsummit']
    df['end'] = df['start'] + 1
    return df.loc[['chrom', 'start', 'end']]
