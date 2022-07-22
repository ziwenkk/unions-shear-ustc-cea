#!/usr/bin/env python3

"""split_sample.py

Split input sample into equi-populated bins.

:Author: Martin Kilbinger <martin.kilbinger@cea.fr>

:Date: 2022

"""

import sys
import os

import numpy as np

from optparse import OptionParser
import matplotlib.pylab as plt

from statsmodels.distributions.empirical_distribution import ECDF

from astropy.io import ascii
from astropy.table import Table 

from unions_wl import sample
from sp_validation import util


def params_default():

    params = {
        'input_path': 'SDSS_SMBH_202206.txt',
        'key_ra': 'ra',
        'key_dec': 'dec',
        'key_z': 'z',
        'key_logM': 'logM',
        'n_split': 2,
        'n_bin_z_hist': 100,
        'output_dir': 'data_mass_sub',
        'output_fname_base': 'SDSS_SMBH_202206',
    }

    short_options = {
        'input_path': '-i',
        'n_split': '-n',
        'output_dir': '-o',
    }

    types = {
        'n_split': 'int',
        'n_bin_z_hist': 'int',
    }

    help_strings = {
        'input_path': 'catalogue input path, default={}',
        'key_ra': 'column name for right ascension, default={}',
        'key_dec': 'column name for declination, default={}',
        'n_spit': 'number of equi-populated bins on output, default={}',
        'n_bin_z_hist': 'number of bins for redshift histogram, default={}',
        'output_dir': 'output directory, default={}',
    }


    return params, short_options, types, help_strings


def parse_options(p_def, short_options, types, help_strings):
    """Parse command line options.

    Parameters
    ----------
    p_def : dict
        default parameter values
    help_strings : dict
        help strings for options

    Returns
    -------
    options: tuple
        Command line options
    """

    usage  = "%prog [OPTIONS]"
    parser = OptionParser(usage=usage)

    for key in p_def:
        if key in help_strings:

            if key in short_options:
                short = short_options[key]
            else:
                short = ''

            if key in types:
                typ = types[key]
            else:
                typ = 'string'

            parser.add_option(
                short,
                f'--{key}',
                dest=key,
                type=typ,
                default=p_def[key],
                help=help_strings[key].format(p_def[key]),
            )

    parser.add_option(
        '-v',
        '--verbose',
        dest='verbose',
        action='store_true',
        help=f'verbose output'
    )

    options, args = parser.parse_args()

    return options


def main(argv=None):

    params, short_options, types, help_strings  = params_default()

    options = parse_options(params, short_options, types, help_strings)

    # Update parameter values
    for key in vars(options):
        params[key] = getattr(options, key)

    # Save calling command
    util.log_command(argv)

    # Open input catalogue
    if params['verbose']:
        print(f'Reading catalogue {params["input_path"]}...')
    dat = ascii.read(f'{path_to_data}/{cat_name}', names=names)

    # To split into more equi-populated bins, compute cumulative distribution function
    cdf = ECDF(dat['logM'])

    # Split into two (check whether we get median from before)
    logM_bounds = sample.y_equi(cdf, params['n_split'])

    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv))
