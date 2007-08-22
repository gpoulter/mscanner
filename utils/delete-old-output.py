#!/usr/bin/env python

from mscanner.configuration import rc
from mscanner.scorefile import readDescriptor

max_keep = 100

if __name__ == "__main__":
    
    # Read outupt descriptors and sort by date
    output = []
    for fpath in rc.web_report_dir.dirs():
        if (fpath / rc.report_index).exists():
            descriptor = readDescriptor(fpath / rc.report_descriptor)
            output.append(descriptor.timestamp, fpath)
    output.sort(reverse=True)
    
    # If more than 100 outpus, delete the old ones
    for timestamp, outpath in output[max_keep:]:
        # Remove files in the 
        for fname in outpath.files():
            fname.remove()
        output.rmdir()