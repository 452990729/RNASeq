#!/usr/bin/env python2


import sys
import re
import os
import ConfigParser
import argparse

BasePath = os.path.split(os.path.realpath(__file__))[0]
config = ConfigParser.ConfigParser()
config.read(BasePath+'/../Bin/config.ini')

#### SOFT
Rscript = config.get('SOFTWARE', 'Rscript')


### SCRIPT
DiffAnalysis = config.get('SCRIPT', 'DiffAnalysis')

#### DATABASE

class ReadList(object):
    def __init__(self, line_in):
        list_split = re.split('\s+', line_in)
        self.Sample = list_split[0]
        self.Group = list_split[1]
        self.Name = '_'.join(list_split[:2])
        self.Line = line_in

def GetSub(compare, dict_group, outpath, raw_path):
    compare_group = re.split(':', compare)
    outdata_dir = os.path.join(outpath, 'VS'.join(compare_group))
    if os.path.exists(outdata_dir):
        os.system('rm -rf {}'.format(outdata_dir))
    os.mkdir(outdata_dir)
    out = open(os.path.join(outdata_dir, 'VS'.join(compare_group)+'.lst'), 'w')
    out_rawdata_dir = os.path.join(outdata_dir, 'rawdata')
    if os.path.exists(out_rawdata_dir):
        os.system('rm -rf {}'.format(out_rawdata_dir))
    os.mkdir(out_rawdata_dir)
    for i in compare_group:
        for ob in dict_group[i]:
            os.system('ln -s {} {}'.format(os.path.join(raw_path, ob.Name),\
                                      os.path.join(out_rawdata_dir, ob.Name)))
            out.write(ob.Line+'\n')
    out.close()
    os.system(' '.join([Rscript, DiffAnalysis, out_rawdata_dir, \
                       os.path.join(outdata_dir, 'VS'.join(compare_group)+'.lst'),\
                       outdata_dir]))

def GetReads(file_in):
    dict_group = {}
    with open(file_in, 'r') as f:
        for line in f:
            ob = ReadList(line.strip())
            if ob.Group not in dict_group:
                dict_group[ob.Group] = [ob, ]
            else:
                dict_group[ob.Group] += [ob, ]
    return dict_group

def main():
    parser = argparse.ArgumentParser(description="Diff Analysis pipeline")
    parser.add_argument('-c', help='the input fasta list', required=True)
    parser.add_argument('-o', help='the abs output path', required=True)
    parser.add_argument('-i', help='the abs input path', required=True)
    parser.add_argument('-g', help='the compare group, A:B,B:C', required=True)
    argv=vars(parser.parse_args())
    dict_group = GetReads(argv['c'])
    for i in re.split(',', argv['g']):
        GetSub(i, dict_group, argv['o'], argv['i'])


if __name__ == '__main__':
    main()




