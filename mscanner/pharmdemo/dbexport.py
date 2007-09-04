"""Exports a database of gene-drug interactions for pharmdemo.stanford.edu

@var schema: Oracle SQL table schema of the pharmdemo database
"""

                                     
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
        self.gdcounts = self.countGeneDrug(gdarticles)

    @staticmethod
    def countGeneDrug(articles):
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
    
    def writeGeneDrugCountsCSV(self, fname):
        """Write PubMed IDs gene-drug associations as CSV
        
        @param fname: Path to output file."""
        import codecs
        f = codecs.open(fname, "wb", "utf-8")
        f.write("PMID,GENE,DRUG\n")
        for (gene,drug),pmids in self.gdcounts.iteritems():
            for pmid in pmids:
                f.write("%s,%s,%s\n" % (pmid, gene, drug))
        f.close()
    
    def exportDatabase(self, con):
        """Export articles to a database connection
    
        @param con: Connection to destination database. Database should 
        be empty before calling this function."""
        cur = con.cursor()
        # Create tables
        cur.executescript(schema)
        # Insert gene-drug associations
        for dg_id, (gd,pmids) in enumerate(self.gdcounts.iteritems()):
            cur.execute('INSERT INTO genedrug(id,gene,drug,numarticles) VALUES (?,?,?,?)',
                        (dg_id, gd[0], gd[1], len(pmids)))
            cur.executemany('INSERT INTO dg_pmids(pmid,text,dg_id) VALUES (?,?,?)',
                            ((str(pmid), None, dg_id) for pmid in pmids))
        # Insert articles
        for art in self.gdarticles:
            genelist = set()
            for genes in art.genedrug.values():
                genelist.update(genes)
            druglist = art.genedrug.keys()
            cur.execute('INSERT INTO cbs(pmid,title,abs,genes,drugs,coes,evid_loc) VALUES (?,?,?,?,?,?,?)',
                        (str(art.pmid), art.title, art.abstract, 
                         " ".join(genelist), " ".join(druglist), None, None))
        con.commit()
    
    def exportText(self, outfile):
        """Export articles to SQL text"""
        
        class TextOutput:
            def __init__(self, outfile):
                import codecs
                self.out = codecs.open(outfile, "wb", "utf-8")
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
        self.exportDatabase(TextOutput(outfile))
    
    def exportSQLite(self, outfile):
        """Export article results to an SQLite database
    
        @param outfile: Path of SQLite database write to."""
        from pysqlite2 import dbapi2 as sqlite
        if outfile.isfile():
            outfile.remove()
        con = sqlite.connect(outfile)
        self.exportDatabase(con)
        con.close()
    
    def exportOracleCon(self, conpath):
        """Export article results to an Oracle database.
    
        @param conpath: user/password@host connection string for the database"""
        import DCOracle2
        con = DCOracle2.connect(conpath)
        self.exportDatabase(cont)
        con.close()

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
