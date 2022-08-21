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
         id_type="library_id",
         index="library_id",
         outfile="celldb.dir/out.load"):
    '''load a table or set of tables into the celldb database'''

    # 1. concatenate separate tables if
    #    a directory has been passed
    #
    if not os.path.isfile(table_path):
        '''concatenate a set of tables'''

        table_file = outfile.replace(".load", ".tsv.gz")

        regex_filename = glob.replace("*",".*/(.*)")

        statement = '''python -m cgatcore.tables
                               --cat %(id_type)s
                               --regex-filename "%(regex_filename)s"
                               %(table_path)s/%(glob)s
                       | gzip -c
                       > %(table_file)s
                    '''

        P.run(statement)

    else:
        table_file = table_path

    # 2. load the table into the database
    #
    idx_stat = " ".join([ "-i " + x + " "
                 for x in index.split(",")])

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
    

    # Add the cell_id index
    # the cell_id index is made on columns "barcode" and "library_id" 
    # some code is taken from https://github.com/cgat-developers/cgat-core/blob/master/cgatcore/csv2db.py
    
    # (i) connect to db

    flavour = csv2db.get_flavour(db_url)
    tablename = csv2db.quote_tablename(table_name,flavour=flavour)
    dbhandle = database.connect(url=db_url)
    
    # (ii) check the table exists and fetch the column names
    
    if not tablename in database.getTables(dbhandle):
    
        raise ValueError("table: %s not found in the database")
    
    # (iii) if "barcode" and "libary_id" columns are present
    # add the cell_id index.
    columns = getColumnNames(dbhandle, tablename)
    
    if "library_id" in columns and "barcode" in columns:

        idx_name = tablename + "_cell_id"

        try:
            statement = "DROP index IF EXISTS " + idx_name
            cc = database.executewait(dbhandle, statement, retries=20)
            cc.close()
            statement = "CREATE index %s ON %s(library_id, barcode)" % (idx_name, tablename)
            cc = database.executewait(dbhandle, statement, retries=20)
            cc.close()
            
        except Exception as ex:
            E.info("adding index %s failed: %s" % (idx_name, ex))
            raise ValueError("Failed to add the idex")
            
        
        
        