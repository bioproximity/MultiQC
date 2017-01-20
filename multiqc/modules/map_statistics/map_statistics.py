#!/usr/bin/env python

""" MultiQC module to parse logs from Skewer """

from __future__ import print_function

import os
from collections import OrderedDict
import logging
import re
from multiqc import config
from multiqc.plots import linegraph
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """ Skewer """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='MapStatistics', anchor='map_statistics',
                                            href="https://?",
                                            info="some test parser")
        self.map_statistics_data = dict()
        self.map_statistics_readlen_dist = dict()

        for f in self.find_log_files(config.sp['map_statistics'], filehandles=True):
            self.parse_map_statistics_log(f)

        if len(self.map_statistics_data) == 0:
            log.debug("Could not find any data in {}".format(config.analysis_dir))
            raise UserWarning

        headers = OrderedDict()
        headers['nr_features'] = {
            'title': 'Number of features',
            'description': '% of reads trimmed',
            'scale': 'RdYlGn-rev',
            'max': 100,
            'min': 0,
            'suffix': '',
        }

        all_hdrs = """
        #slice	RT_begin	RT_end	number_of_features	tic	int_mean
        int_stddev	int_min	int_max	int_median	int_lowerq
        int_upperq	mz_mean	mz_stddev	mz_min	mz_max	mz_median
        mz_lowerq	mz_upperq	width_mean	width_stddev
        width_min	width_max	width_median	width_lowerq
        width_upperq	qual_mean	qual_stddev	qual_min	qual_max
        qual_median	qual_lowerq	qual_upperq	rt_qual_mean
        rt_qual_stddev	rt_qual_min	rt_qual_max	rt_qual_median
        rt_qual_lowerq	rt_qual_upperq	mz_qual_mean	mz_qual_stddev
        mz_qual_min	mz_qual_max	mz_qual_median	mz_qual_lowerq
        mz_qual_upperq
        """
        for h in all_hdrs.split():
            headers[h] = {
                'title': h,
                'description': '?',
                'format': '{.3e}',
            }



        self.general_stats_addcols(self.map_statistics_data, headers)

        # Write parsed report data to a file
        self.write_data_file(self.map_statistics_data, 'multiqc_map_statistics')
        self.write_data_file(self.map_statistics_readlen_dist, 'multiqc_map_statistics_readlen_dist')

        # set the value 0 for every x where a given sample doens't have a value
        all_x_values = []
        for s_name in self.map_statistics_readlen_dist:
            for xval in self.map_statistics_readlen_dist[s_name]:
                all_x_values.append(xval)

        for s_name in self.map_statistics_readlen_dist:
            for xval in all_x_values:
                if not xval in self.map_statistics_readlen_dist[s_name]:
                    self.map_statistics_readlen_dist[s_name][xval] = 0.0

            # After adding new elements, the ordereddict needs to be re-sorted
            items = self.map_statistics_readlen_dist[s_name]
            self.map_statistics_readlen_dist[s_name] = OrderedDict(sorted(items.items(), key=lambda x: int(x[0])))

        log.info("Found {} reports".format(len(self.map_statistics_data)))

    def parse_map_statistics_log(self, f):
        """ Go through log file looking for map_statistics output """
        fh = f['f']
        regexes = {
            'nr_features': "Number of features: (\d+)",
        }

        data = dict()
        for k, v in regexes.items():
            data[k] = 0
        data['nr_features'] = None

        while True:
            l=fh.readline()
            if l=="": break
            for k, r in regexes.items():
                match = re.search(r, l)
                if match:
                    data[k] = match.group(1)
            if l.startswith("-- Summary Statistics --"):
                fh.readline()
                hdrs=fh.readline().rstrip("\n").split("\t")
                cols=[ float(a) for a in fh.readline().rstrip("\n").split("\t") ]
                d = dict(zip(hdrs, cols))
                data.update(d)


        if data['nr_features'] is not None:
            s_name = self.clean_s_name(f['s_name'], f['root'])
            if s_name in self.map_statistics_data:
                log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
            self.map_statistics_data[s_name] = {}
            self.map_statistics_data[s_name] = data.copy()
