#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import datetime

'''
Fastq_file_prefix
ReadGroup_sample_RGSM
ReadGroup_id_RGID
ReadGroup_library_RGLB
ReadGroup_platform_RGPL
ReadGroup_platform_unit_RGPU
ReadGroup_SeqCentre_RGCN
ReadGroup_Desc_RGDS
ReadGroup_runDate_RGDT
QUAL
Pipeline
PE
bed_list
bed_type
email_bioinf
Fastq_dir
BAM_dir
pipeline_dir
sleep
'''

class Patients(list):
    def __init__(self, fh, key=None):
        # the sorting key provides some flexibility in ordering the features
        self.key = key or (lambda x: (x.chrom, x.chromStart, x.chromEnd))
        for line in fh:
            if line[0] in ("#", "/") or line.startswith('Fastq_file_prefix'):
                continue
            self.append(BEDline(line))

        self.seqids = sorted(set(b.chrom for b in self))
        self.sort(key=self.key)
        return

    def get_order(self):
        return dict((f.name, (i, f)) for (i, f) in enumerate(self))

    def get_simple_bed(self):
        return [(b.chrom, i) for (i, b) in enumerate(self)]

    def getDict(self, keyattribute='name'):
        selfdict = {}
        for b in self:
            selfdict[getattr(b,keyattribute)] = b
        return selfdict



class Patient(list):

    __slots__ = ("Fastq_file_prefix","ReadGroup_sample_RGSM",
        "ReadGroup_id_RGID","ReadGroup_library_RGLB",
        "ReadGroup_platform_RGPL","ReadGroup_platform_unit_RGPU",
        "ReadGroup_SeqCentre_RGCN","ReadGroup_Desc_RGDS",
        "ReadGroup_runDate_RGDT","QUAL","Pipeline","PE","bed_list",
        "bed_type","email_bioinf","Fastq_dir","BAM_dir","pipeline_dir",
        "sleep")


    def __init__(self, sline):
        args = sline.strip().split("\t")
        try:
            assert len(args) == 19
        except:
            print >> sys.stderr, '##', sline, '##'
            raise Exception('Patient data parsing error')
        self.meta = {}
        # read fields
        for i in range(len(args)):
            setattr(self, Patient.__slots__[i], args[i])
    return

    def __str__(self):
        fields = []
        for i in range(len(Patient.__slots__)):
            try:
                fields.append(getattr(self, Patient.__slots__[i]))
            except:
                break
        return "\t".join(map(str, fields))

    def __getitem__(self, key):
        return getattr(self, key)

    def __lt__(self, other):  # sort by time year-month-day
        try:
            selfDate = tuple(map(int, self.ReadGroup_runDate_RGDT.split("-")[:3]))
            otherDate = tuple(map(int, other.ReadGroup_runDate_RGDT.split("-")[:3])) 
        except:
            logging.warning("date sorting error: not YYYY-MM-DD")
            return True
        return selfDate < otherDate


if __name__ == "__main__":
    patients = Patients(sys.stdin)
    for patient in patient:
        print patient