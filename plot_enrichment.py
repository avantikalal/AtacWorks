#
# Copyright (c) 2019, NVIDIA CORPORATION.  All rights reserved.
#
# NVIDIA CORPORATION and its licensors retain all intellectual property
# and proprietary rights in and to this software, related documentation
# and any modifications thereto.  Any use, reproduction, disclosure or
# distribution of this software and related documentation without an express
# license agreement from NVIDIA CORPORATION is strictly prohibited.
#

"""Create an enrichment plot from a bigWig file.

"""

# import requirements

import numpy as np
import matplotlib.pyplot as plt
import argparse
import logging
from claragenomics.io.bigwigio import extract_bigwig_intervals
from claragenomics.io.bedio import read_intervals, read_summits, read_sizes

# Set up logging
log_formatter = logging.Formatter(
    '%(levelname)s:%(asctime)s:%(name)s] %(message)s')
_logger = logging.getLogger('AtacWorks-plot-tss')
_handler = logging.StreamHandler()
_handler.setLevel(logging.INFO)
_handler.setFormatter(log_formatter)
_logger.setLevel(logging.INFO)
_logger.addHandler(_handler)


def reverse_negative_strand(x, strand):
    """Reverse the direction of a signal on the negative strand.

    Args:
        x: 2-D numpy array containing signal in rows
        strand: list of strand ('+' or '-') for each row
    Return:
        x_reversed: numpy array,
        same as x with signal reversed when strand == "-".

    """
    num_signals = x.shape[0]
    assert num_signals == len(strand)
    x_reversed = [
        np.flip(x[i, :]) if strand[i] == '-' else x[i, :] for i in range(
            num_signals)]
    x_reversed = np.stack(x_reversed)
    assert x_reversed.shape == x.shape
    return x_reversed


def parse_args():
    """Parse command line arguments.

    Return:
        parsed argument object.

    """
    parser = argparse.ArgumentParser(
        description='Script to plot enrichment in bigWig')
    parser.add_argument('--bw', type=str, required=True,
                        help='Path to bigWig file')
    parser.add_argument('--bed', type=str,
                        help='Path to BED file with genomic sites')
    parser.add_argument('--summits', type=str,
                        help='Path to narrowPeak file with summits')
    parser.add_argument('--sizes', type=str,
                        help='Path to chromosome sizes file.')
    parser.add_argument('--prefix', type=str, help='Output file prefix',
                        required=True)
    parser.add_argument('--extend', type=int,
                        help='Extension on either side of genomic site',
                        default=2000)
    parser.add_argument('--stranded', action='store_true',
                        help='Intervals are stranded')
    parser.add_argument('--norm', type=str,
                        help='range of positions by which to normalize \
                        signal, start:end', default='1:200')
    args = parser.parse_args()
    return args


args = parse_args()

# Read genomic sites or summits
if args.bed is not None:
    _logger.info('Reading BED file')
    sites = read_intervals(args.bed, stranded=args.stranded)
elif args.summits is not None:
    _logger.info('Reading narrowPeak file')
    sites = read_summits(args.summits)
_logger.info('Read ' + str(len(sites)) + ' genomic sites.')

# Extend
_logger.info('Extending sites by ' + str(args.extend) + ' on either side.')
sites['start'] = sites['start'] - args.extend
sites['end'] = sites['end'] + args.extend

# Read sizes
if args.sizes is not None:
    _logger.info('Reading chromosome sizes file')
    sizes = read_sizes(args.sizes)

    # Filter chromosomes
    sites_out = sites[~sites['chrom'].isin(sizes['chrom'])]
    _logger.info("Discarding " + str(
        len(sites_out)) + " sites in chromosomes outside the sizes file")
    sites = sites[sites['chrom'].isin(sizes['chrom'])]

# Read data
_logger.info('Extracting signal around ' + str(len(sites)) + ' genomic sites')
signal = extract_bigwig_intervals(sites, args.bw)

# Flip values for genes on negative strand
if args.stranded:
    _logger.info('Reversing signals on negative strand')
    signal = reverse_negative_strand(signal, sites['strand'])

# Sum
_logger.info('Adding signal over all sites')
signal = np.sum(signal, axis=0)[:-1]

# Normalize
_logger.info('Normalizing signal relative to positions ' + args.norm)
norm_start, norm_end = args.norm.split(':')
norm_factor = signal[int(norm_start):int(norm_end)].mean()
signal = signal / norm_factor

# Plot
_logger.info('Saving plot to ' + args.prefix + '.pdf')
plt.plot(range(-1 * args.extend, args.extend), signal)
plt.savefig(args.prefix + '.pdf')
