#!/usr/bin/env python
#
# Copyright (c) 2019, NVIDIA CORPORATION.  All rights reserved.
#
# NVIDIA CORPORATION and its licensors retain all intellectual property
# and proprietary rights in and to this software, related documentation
# and any modifications thereto.  Any use, reproduction, disclosure or
# distribution of this software and related documentation without an express
# license agreement from NVIDIA CORPORATION is strictly prohibited.
#

"""Returns quality metrics for a coverage track in bigWig format

"""

# Import requirements
import argparse
import numpy as np
import logging
from claragenomics.io.bigwigio import extract_bigwig_intervals
from claragenomics.io.bedio import read_sizes, read_summits, read_intervals

# Set up logging
log_formatter = logging.Formatter(
    '%(levelname)s:%(asctime)s:%(name)s] %(message)s')
_logger = logging.getLogger('AtacWorks-qual-metrics')
_handler = logging.StreamHandler()
_handler.setLevel(logging.INFO)
_handler.setFormatter(log_formatter)
_logger.setLevel(logging.INFO)
_logger.addHandler(_handler)


def parse_args():
    """Parse command line arguments.

    Return:
        parsed argument object.

    """
    parser = argparse.ArgumentParser(
        description='Calculate quality metrics from bigWig file and peaks')
    parser.add_argument('--bw', type=str, required=True,
                        help='Path to bigWig file')
    parser.add_argument('--sizes', type=str, required=True,
                        help='Path to chromosome sizes file.')
    parser.add_argument('--peaks', type=str, required=True,
                        help='Path to narrowPeak file with summits')
    args = parser.parse_args()
    return args


args = parse_args()

# Read sizes
_logger.info('Loading chromosome sizes')
sizes = read_sizes(args.sizes, as_intervals=True)

# Read peaks
_logger.info('Loading peaks')
peaks = read_intervals(args.peaks, skip=1)

# Read summits
_logger.info('Loading summits')
summits = read_summits(args.peaks)

# Filter peaks
peaks_drop = peaks[~peaks['chrom'].isin(sizes['chrom'])]
peaks_keep = peaks['chrom'].isin(sizes['chrom'])
_logger.info('Discarding ' + str(len(
    peaks_drop)) + ' peaks outside sizes file.')
peaks = peaks[peaks_keep]
peaks.reset_index(drop=True, inplace=True)

# Filter summits
summits = summits[peaks_keep]
summits.reset_index(drop=True, inplace=True)

# Number of peaks
_logger.info('Counting peaks')
num_peaks = len(summits)
print("Number of peaks: {}".format(num_peaks))

# Extract overall signal
_logger.info('Extracting signal values for all positions')
signal = extract_bigwig_intervals(sizes, args.bw, stack=False)
signal = np.concatenate(signal)

# Calculate mean signal
_logger.info('Calculating average signal')
sum_signal = np.sum(signal)
mean_signal = sum_signal / len(signal)

# FC at summits
_logger.info('Extracting signal at summits')
signal_at_summits = extract_bigwig_intervals(summits, args.bw)
signal_at_summits = np.concatenate(signal_at_summits)
fc_at_summits = signal_at_summits / mean_signal
assert(len(fc_at_summits) == num_peaks)

# Number of peaks with fold enrichment over threshold
_logger.info('Counting peaks with fold enrichment over threshold')
num_peaks_10 = sum(fc_at_summits >= 10)
print("Number of peaks with FC>=10 over global average: {}".format(
    num_peaks_10))
num_peaks_20 = sum(fc_at_summits >= 20)
print("Number of peaks with FC>=20 over global average: {}".format(
    num_peaks_20))

# FSIP
_logger.info('Calculating FSIP')
signal_in_peaks = extract_bigwig_intervals(peaks, args.bw, stack=False)
signal_in_peaks = np.concatenate(signal_in_peaks)
fsip = np.sum(signal_in_peaks) / sum_signal
print("FSIP: {}".format(fsip))
