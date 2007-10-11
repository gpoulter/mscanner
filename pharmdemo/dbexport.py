"""Exports a database of gene-drug interactions for pharmdemo.stanford.edu

@var schema: Oracle SQL table schema of the pharmdemo database
"""

from __future__ import with_statement
from contextlib import closing


                                     
__author__ = "Graham Poulter"                                        
__license__ = """This program is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your option)
any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>."""


class GeneDrugExport:
    """Exports database of gene-drug interactions 
    
    @ivar gdarticles: List of articles (with additional genedrug members)
    @ivar gdcounts: Mapping of (gene,drug) pairs to list of PubMed IDs
    """
    
    def __init__(self, gdarticles):
        self.gdarticles = gdarticles
        self.gdcounts = self.count_genedrug(gdarticles)


    @staticmethod
    def count_genedrug(articles):
        """Return global gene-drug associations given per-article associations
        
        @param articles: List of Article objects. Requires pmid, title,
        abstract and genedrug fields. with genedrug field
        
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


    def write_genedrug_csv(self, fname):
        """Write PubMed IDs gene-drug associations as CSV
        
        @param fname: Path to output file."""
        import codecs
        with codecs.open(fname, "wb", "utf-8") as f:
            f.write("PMID,GENE,DRUG\n")
            for (gene,drug),pmids in self.gdcounts.iteritems():
                for pmid in pmids:
                    f.write("%s,%s,%s\n" % (pmid, gene, drug))


    def export_dbapi(self, con):
        """Export articles to a database connection
        
        @param con: Connection to destination database. Database should 
        be empty before calling this function."""
        cur = con.cursor()
        cur.executescript(schema) # Create tables 
        # Insert gene-drug associations
        for dg_id, (gd,pmids) in enumerate(self.gdcounts.iteritems()):
            # INSERT INTO genedrug(id,gene,drug,numarticles) 
            cur.execute("INSERT INTO genedrug VALUES (?,?,?,?)",
                        (dg_id, gd[0], gd[1], len(pmids)))
            # INSERT INTO dg_pmids(pmid,text,dg_id) 
            cur.executemany("INSERT INTO dg_pmids VALUES (?,?,?)",
                            ((str(pmid), None, dg_id) for pmid in pmids))
        # Insert citations
        for art in self.gdarticles:
            genelist = set()
            for genes in art.genedrug.values():
                genelist.update(genes)
            druglist = art.genedrug.keys()
            # INSERT INTO cbs(pmid,title,abs,genes,drugs,coes,evid_loc) 
            # The '\n' is to put the abstract on its own line when
            # exporting as SQL text, as sometimes the whole query goes
            # over the Oracle sqplus 2499 character per line limit.
            abstract = art.abstract[:2450] if art.abstract else None
            cur.execute("INSERT INTO cbs VALUES (?,?,\n?,\n?,?,?,?)",
                        (str(art.pmid), art.title, abstract, 
                         " ".join(genelist), " ".join(druglist), None, None))
        con.commit()


    def export_sqlfile(self, outfile):
        """Export to text file"""
        class TextOutput:
            """'Connection' that just writes SQL to file
            Oracle single quotes go around string data.            
            """
            def __init__(self, outfile):
                import codecs
                self.out = codecs.open(outfile, "wb", "utf-8")
            def commit(self):
                self.out.flush()
            def close(self):
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
                    elif isinstance(s, basestring):
                        sub[i] = "'" + s.replace("'","''") + "'"
                    else:
                        sub[i] = str(s)
                self.out.write(sql % tuple(sub))
            def executemany(self, sql, subs):
                for sub in subs:
                    self.execute(sql, sub)
        with closing(TextOutput(outfile)) as con:
            self.export_dbapi(con)


    def export_sqlite(self, outfile):
        """Export to an SQLite database
        
        @param outfile: Path of SQLite database write to."""
        from pysqlite2 import dbapi2 as sqlite
        if outfile.isfile():
            outfile.remove()
        with closing(sqlite.connect(outfile)) as con:
            self.export_dbapi(con)


    def export_dcoracle2(self, conpath):
        """Export to an Oracle database.
        
        @param conpath: user/password@host connection string for the database"""
        import DCOracle2
        with closing(DCOracle2.connect(conpath)) as con:
            self.export_dbapi(cont)



schema = """
-- suppress sql in result set
set echo off
-- eliminate row count message
set feedback off
-- disable ampersand substitution
set define off
-- suppress headings and page breaks
set pagesize 0
-- suppress all terminal echo
set term off

-- Drop tables in reverse order
drop table dg_pmids;
drop table genedrug;
drop table cbs;

-- Articles, categories of evidence, and genes/drugs
create table cbs (
  pmid varchar2(20) primary key,
  title varchar2(2000),
  abs clob,
  genes varchar2(2000),
  drugs varchar2(2000),
  coes varchar2(50),
  evid_loc varchar2(50)
);

-- Gene-drug pairs, and number of articles containing both
create table genedrug (
  id number(10,0) primary key,
  gene varchar2(100),
  drug varchar2(100),
  numarticles number(10,0)
);

-- PMIDs associated with the gene-drug pairs
create table dg_pmids (
  pmid varchar2(15),
  text clob,
  dg_id integer,
  foreign key (dg_id) references genedrug(id)
);
"""
"""Database schema prepended to the exports"""


schema_unused = """
drop table dg_pmid_coes;

-- Category of evidence associated with each article
create table dg_pmid_coes (
  pmid varchar2(15) primary key,
  coe varchar2(3)
);
"""
"""Unused schema - since we're not doing category of evidence"""