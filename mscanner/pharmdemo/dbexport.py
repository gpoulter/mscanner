"""Export a database compatible with pharmdemo.stanford.edu

countGeneDrug() -- Turn per-article associations into global associations
writeGeneDrugCountsCSV() -- Record global gene-drug associations for checking
exportDatabase() -- Export articles to a database connection
exportText() -- Export articles to SQL text
exportSQLite() -- Export articles to an SQLite database
schema -- The schema of the pharmdemo database

                                   

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

http://www.gnu.org/copyleft/gpl.html
"""

import logging as log
import csv
import codecs

class ucsvwriter:
    """A CSV writer taking unicode objects, writing to a stream in utf-8"""
    def __init__(self, f, dialect=csv.excel, **kwds):
        self.writer = codecs.getwriter(f, dialect=dialect, **kwds)
        self.dialect = dialect
    def writerow(self, row):
        self.writer.writerow([s.encode("utf-8") for s in row])
    def writerows(self, rows):
        for row in rows:
            self.writerow(row)

def ucsvreader(f, dialect=csv.excel, **kwds):
    """A CSV reader for a stream in utf-8, returning unicode objects"""
    for line in csv.reader(f, dialect=dialect, **kwds):
        return unicode(line, "utf-8")

def countGeneDrug(articles):
    """Return global gene-drug associations given per-article associations

    @param articles: Iterable of Article objects with genedrug field

    @return: Mapping from (gene,drug) pair to list of PMIDs
    """
    # Count gene-drug co-occurrences from articles
    gdcounter = {}
    for art in articles:
        for drug, genes in art.genedrug.iteritems():
            for gene in genes:
                if (gene,drug) not in gdcounter:
                    gdcounter[(gene,drug)] = []
                gdcounter[(gene,drug)].append(art.pmid)
    return gdcounter

def writeGeneDrugCountsCSV(fname, gdcounts):
    """Write association of PubMed IDs to gene-drug occurrences to a CSV file"""
    outf = ucsvwriter(file("/tmp/output.csv", "wb"))
    outf.writerow("PMID,GENE,DRUG\n")
    for (gene,drug),pmids in gdcount.iteritems():
        for pmid in pmids:
            outf.writerow("%s,%s,%s\n" % (pmid, gene, drug))

def exportDatabase(con, articles):
    """Export articles to database

    @param con: Connection to destination database. Database should be
    empty.

    @param articles: Articles to export.  Requires pmid, title,
    abstract and genedrug fields.
    """
    cur = con.cursor()
    # Create tables
    cur.executescript(schema)
    # Insert gene-drug associations
    gdcounter = countGeneDrug(articles)
    for dg_id, (gd,pmids) in enumerate(gdcounter.iteritems()):
        cur.execute('INSERT INTO genedrug(id,gene,drug,numarticles) VALUES (?,?,?,?)',
                    (dg_id, gd[0], gd[1], len(pmids)))
        cur.executemany('INSERT INTO dg_pmids(pmid,text,dg_id) VALUES (?,?,?)',
                        ((str(pmid), None, dg_id) for pmid in pmids))
    # Insert articles
    for art in articles:
        genelist = set()
        for genes in art.genedrug.values():
            genelist.update(genes)
        druglist = art.genedrug.keys()
        cur.execute('INSERT INTO cbs(pmid,title,abs,genes,drugs,coes,evid_loc) VALUES (?,?,?,?,?,?,?)',
                    (str(art.pmid), art.title, art.abstract, " ".join(genelist), " ".join(druglist), None, None))
    con.commit()

def exportText(outfile, articles):
    class TextOutput:
        def __init__(self, outfile):
            self.out = file(outfile, "w")
        def __del__(self):
            self.commit()
        def commit(self):
            self.out.close()
        def cursor(self):
            return self
        def executescript(self, sql):
            self.out.write(sql)
        def execute(self, sql, sub):
            sql = sql.replace("?", "%s") + ";\n"
            sub = list(sub)
            for i,s in enumerate(sub):
                if s is None:
                    sub[i] = "NULL"
                elif isinstance(s, str):
                    sub[i] = '"' + str(s) + '"'
            self.executescript(sql % tuple(sub))
        def executemany(self, sql, subs):
            for sub in subs:
                self.execute(sql, sub)
    exportDatabase(TextOutput(outfile), articles)

def exportSQLite(outfile, articles):
    """Export article results to SQLite database file

    @param outfile: Path of SQLite database write to.

    @param articles: List of Article objects to export.  Requires
    pmid, title, abstract and genedrug fields.
    """
    from pysqlite2 import dbapi2 as sqlite
    log.info("Exporting pharmdemo database to %s", outfile.name)
    if outfile.isfile():
        outfile.remove()
    con = sqlite.connect(outfile)
    exportDatabase(con, articles)
    con.close()

def exportOracleCon(conpath, articles):
    """Export article results to an Oracle database.

    @param conpath: user/password@host for the database
    """
    import DCOracle2
    con = DCOracle2.connect(conpath)
    exportDatabase(cont, articles)
    con.close()

exportDefault = exportText

schema="""
-- lists articles, providing categories of evidence (COE) information
-- as well as unprocessed lists of genes and drugs present in the abstract

create table cbs (
  pmid varchar2(20) primary key,
  title varchar2(2000),
  abs clob,
  genes varchar2(2000),
  drugs varchar2(2000),
  coes varchar2(50),
  evid_loc varchar2(50)
);

-- contains gene-drug pairs, and the number of articles containing
-- both in their abstract.

create table genedrug (
  id number(10,0) primary key, 
  gene varchar2(100), 
  drug varchar2(100),
  numarticles number(10,0) 
);

-- lists PMIDs associated with the gene-drug pairs from the
-- previous table.

create table dg_pmids (
  pmid varchar2(15),
  text clob,
  dg_id integer,
  foreign key (dg_id) references genedrug(id)
);

-- lists category of evidence associated with each article

create table dg_pmid_coes (
  pmid varchar2(15),
  coe varchar2(3)
);
"""