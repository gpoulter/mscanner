#!env python

"""Parses a drug table

Usage: parse_drugs.py <file>

The specified file should map PharmGKB Accession ID's to drug names.  They are
used to scan abstracts for occurrences of drugs.  One drug is
described per line, and each line is to be formatted as follows:

The drug's PharmGKB Accession ID
[tab]
The drug's main name in PharmGKB (might be comma-separated names)
[tab]
A list of the drug's generic names (if any), each terminated by '|'
[tab]
A list of the drug's trade names (if any), each terminated by '|'

The output is 'drugs.cpickle', to be placed in the working directory
of pgmedline.

author: Graham Poulter
                                  
"""

import re, os.path

def parseDrugs(text):
    """Return a drug dictionary from specially formatted table

    @see: Format of the drug table is described in README.html

    @raise ValueError: On incorrect formatting.

    @type text: C{string}
    @param text: Drug table to load.

    @rtype: C{{string:[string]}}
    @return: Mapping from PharmGKB ID to aliases for the drug.  The
    first alias is the preferred drug name.
    """
    drugs = {}
    lines = text.split('\n')
    mainsplit = re.compile(r"\A(\w+)\t([^\t]+)\t([^\t]*)\t(.*)\Z")
    subsplit = re.compile(r"\|")
    for line in lines:
        if line=="":
            continue
        # Check tab count in lines
        if line.count('\t') != 3:
            raise ValueError, "Drug lines must have three tabs (see README). Bad line was '"+line+"'"
        # Check splitting line into components
        (PKID,main,trade,generic) = mainsplit.match(line).groups()
        if trade != "" and trade[-1:] != "|":
            raise ValueError, "Drugs must terminate in '|' (see README). Bad line was '"+line+"'"
        if generic != "" and generic[-1:] != "|":
            raise ValueError, "Drugs must terminate in '|' (see README). Bad line was '"+line+"'"
        # Get names from each component and slurp into a list
        trades = subsplit.split(trade[:-1])
        generics = subsplit.split(generic[:-1])
        drugs[PKID] = [ name for name in [main]+trades+generics if name != "" ]
    return drugs

def testParseDrugs():
    drugs_table = """PA10000\t17 beta-estradiol\t\t
PA10007\talbumin human\t\tAlbuminar-25|Albuminar-5|Albutein 25%|Albutein 5%|Buminate 25%|Buminate 5%|Plasbumin-25|Plasbumin-5|
PA10009\talefacept\t\tAmevive|
PA1001\t1-methyl-4-phenylpyridinium (MPP+)\t\t
PA10010\talemtuzumab\t\tCampath|
PA10011\talfacalcidol\t\tOne-Alpha|
PA10012\talteplase, recombinant\t\tActivase|Activase rt-PA|Cathflo Activase|"""
    fulldrugs_correct = {
        'PA10000': ['17 beta-estradiol'],
        'PA10007': ['albumin human',
                    'Albuminar-25',
                    'Albuminar-5',
                    'Albutein 25%',
                    'Albutein 5%',
                    'Buminate 25%',
                    'Buminate 5%',
                    'Plasbumin-25',
                    'Plasbumin-5'],
        'PA10009': ['alefacept', 'Amevive'],
        'PA1001': ['1-methyl-4-phenylpyridinium (MPP+)'],
        'PA10010': ['alemtuzumab', 'Campath'],
        'PA10011': ['alfacalcidol', 'One-Alpha'],
        'PA10012': ['alteplase, recombinant',
                    'Activase',
                    'Activase rt-PA',
                    'Cathflo Activase']}
    result = parseDrugs( drugs_table )
    assert result == fulldrugs_correct
    
if __name__=="__main__":
    testParseDrugs()
    import sys
    if(len(sys.argv)!=2 or (not os.path.exists(sys.argv[1]))):
        print __doc__
        sys.exit(1)
    drugfile = sys.argv[1]
    drugs = parseDrugs(file(drugfile,'r').read())
    from cPickle import dump
    dump(drugs,file('drugs.cpickle','wb'),protocol=2)
