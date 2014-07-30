#!/usr/bin/env python

import sys
from subprocess import PIPE, call
import logging
from collections import deque

def getColumns(inFile, delim="\t", header=True):
    """
    Get columns of data from inFile. The order of the rows is respected

    :param inFile: column file separated by delim
    :param header: if True the first line will be considered a header line
    :returns: a tuple of 2 dicts (cols, indexToName). cols dict has keys that
    are headings in the inFile, and values are a list of all the entries in that
    column. indexToName dict maps column index to names that are used as keys in
    the cols dict. The names are the same as the headings used in inFile. If
    header is False, then column indices (starting from 0) are used for the
    heading names (i.e. the keys in the cols dict)
    """
    cols = {}
    indexToName = {}
    for lineNum, line in enumerate(inFile):
        if lineNum == 0:
            headings = line.split(delim)
            i = 0
            for heading in headings:
                heading = heading.strip()
                if header:
                    cols[heading] = []
                    indexToName[i] = heading
                else:
                    # in this case the heading is actually just a cell
                    cols[i] = [heading]
                    indexToName[i] = i
                i += 1
        else:
            cells = line.split(delim)
            i = 0
            for cell in cells:
                cell = cell.strip()
                cols[indexToName[i]] += [cell]
                i += 1

    return cols, indexToName

# subprocess with logging
def syscall(*exe):
    try:
        retcode = call(list(exe), shell=True)
        if retcode < 0:
            logging.error("Child was terminated by signal " + str(retcode))
        else:
            print retcode
            logging.debug("Child returned " + str(retcode))
    except OSError as e:
        logging.critical("Execution of "+' '.join(list(exe))+" failed")
        raise e
    return

# call something on SGE
def sgecall(execstring):
    raise Exception("NOT IMPLEMENTED")
    return


# temporary solution until we port everything to drmaa-ruffus
class SGE(object):
    def __init__(self, out=None, err=None, queue=None, email=None, sendemail='a'):
        self._stdout = out
        self._stderr = err
        self._queue = queue
        self._email = email
        self._send = sendemail  # when to send mail (beasn)
        self._resources = []
        self._depends = deque()
        self.optional = []
        return

    def setResources(self, multithread=None, vmem=8):
        self.resources = []
        if multithread:
            self.resources.append('-pe multi_thread %d' % multithread)
        if vmem:
            self.resources.append('-l h_vmem %dG' % vmem)
        return

    def qsub(self, exe, name):
        # global parameters
        command = [ 'qsub',
           '-o '+self._stdout,
           '-e '+self._stderr,
           '-q '+self._queue,
           '-M '+self._email,
           '-m '+self._send ]
        ## add requested resources
        if self._resources:
            command += self._resources
        ## job dependencies
        try:
            command.append('-hold_jid '+ self._depends.popleft())
        except IndexError:
            pass
        # JOB SPECIFIC STUFF
        ## job name
        command.append('-N '+name)
        ## optional
        if self.optional:
            command += self.optional



