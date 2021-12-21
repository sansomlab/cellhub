import os
from cgatcore import pipeline as P

def load(table_name,
         table_path,
         db_url="sqlite:///./celldb.dir/csvdb",
         glob="*.tsv.gz",
         id_type="library_id",
         index="library_id",
         outfile="celldb.dir/out.load",
         mem="8G"):
    '''load a table or set of tables into the celldb database'''

    job_memory = mem

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
