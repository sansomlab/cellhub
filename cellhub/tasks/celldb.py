"""
celldb.py
=========

Helper functions for pipeline_celldb.py

Code
====

"""

import os
from cgatcore import pipeline as P
from cgatcore import database as database
from cgatcore import csv2db as csv2db

def getColumnNames(dbhandle, tablename):

    statement = "PRAGMA table_info(%s);" % tablename
    cc = database.executewait(dbhandle, statement, retries=20)
    columns = [x[1] for x in cc.fetchall()]
    cc.close()
    
    return columns


def load(table_name,
         table_path,
         db_url="sqlite:///./celldb.dir/csvdb",
         glob="*.tsv.gz",
         index=None,
         outfile="celldb.dir/out.load"):
    '''load a table or set of tables into the celldb database'''


    job_memory = "16G"
    job_threads = 1

    # 1. concatenate separate tables if
    #    a directory has been passed
    #
    if not os.path.isfile(table_path):
        '''concatenate a set of tables'''

        table_file = outfile.replace(".load", ".tsv.gz")

        regex_filename = glob.replace("*",".*/(.*)")

        statement = '''python -m cgatcore.tables
                               --regex-filename "%(regex_filename)s"
                               --cat "filename"
                               %(table_path)s/%(glob)s
                       | gzip -c
                       > %(table_file)s
                    '''

        P.run(statement)

    else:
        table_file = table_path

    # 2. load the table into the database
    #
    if index is not None:
        idx_stat = " ".join([ "-i " + x + " "
                 for x in index.split(",")])
    else:
        idx_stat = ""

    # table headers are sanitised to replace
    # spaces and hypens with underscores

    if table_file.endswith(".gz"):
        cat = "zcat"
    else:
        cat = "cat"

    statement='''%(cat)s %(table_file)s
                 | grep -v ^#
                 | sed '1!b;s/[ -]/_/g'
                 | python -m cgatcore.csv2db
                          --retry
                          --database-url=%(db_url)s
                          --table=%(table_name)s
                          %(idx_stat)s
                 > %(outfile)s
              '''

    P.run(statement)
    
    # Create indexes on expected columns
    
    # Connect to the db
    flavour = csv2db.get_flavour(db_url)
    tablename = csv2db.quote_tablename(table_name,flavour=flavour)
    dbhandle = database.connect(url=db_url)
    
    # Check the table exists and fetch the column names
    if not tablename in database.getTables(dbhandle):
    
        raise ValueError("table: %s not found in the database")
    

    columns = getColumnNames(dbhandle, tablename)
        
    def _addIndex(dbhandle, tablename, index_name, index_columns):
    
        try:
            statement = "DROP index IF EXISTS " + index_name
            cc = database.executewait(dbhandle, statement, retries=20)
            cc.close()
            statement = "CREATE index %s ON %s(%s)" % (index_name, 
                                                       tablename, 
                                                       index_columns)
            
            cc = database.executewait(dbhandle, statement, retries=20)
            cc.close()
            
        except Exception as ex:
            print("adding index %s failed: %s" % (index_name, ex))
            raise ValueError("Failed to add the index")
        
    # Add cell_id index
    if "library_id" in columns and "barcode" in columns:

        idx_name = tablename + "_cell_id"
        _addIndex(dbhandle, tablename, idx_name, "library_id, barcode")

    # Add library_id index
    if "library_id" in columns:

        idx_name = tablename + "_library_id"
        _addIndex(dbhandle, tablename, idx_name, "library_id")
    
    # Add sample_id index 
    if "sample_id" in columns:

        idx_name = tablename + "_sample_id"
        _addIndex(dbhandle, tablename, idx_name, "sample_id")
