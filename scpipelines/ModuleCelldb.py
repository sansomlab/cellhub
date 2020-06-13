connection=apsw.Connection("csvdb")
cursor=connection.cursor()

def getfiledata(path):
    columns=None
    data=[]
    counter=1
    for p in path:
        with open(p, "r") as infile:
            for line in infile:
                counter+=1

                if columns is None:
                    columns= line

                else:
                    data.append(line.replace("\t", ",").strip().split(','))
    return columns, data



# This gets registered with the Connection
class Source:
    def Create(self, db, modulename, dbname, tablename, *args):
        columns,data=getfiledata([x for x in args])# eval strips off layer of quotes
        columns = ['%s' % (x,) for x in columns.split()]
        schema="create table foo("+ ','.join(["'%s'" %  x for x in columns]) +")"

        return schema,Table(columns,data)
    Connect=Create

# Represents a table
class Table:
    def __init__(self, columns, data):
        self.columns=columns
        self.data=data

    def BestIndex(self, *args):
        return None

    def Open(self):
        return Cursor(self)

    def Disconnect(self):
        pass

    Destroy=Disconnect

# Represents a cursor
class Cursor:
    def __init__(self, table):
        self.table=table

    def Filter(self, *args):
        self.pos=0

    def Eof(self):
        return self.pos>=len(self.table.data)

    def Rowid(self):
        return self.table.data[self.pos][0]

    def Column(self, col):

        return self.table.data[self.pos][col]

    def Next(self):
        self.pos+=1

    def Close(self):
        pass


connection.createmodule("tsv", Source())
sysdirs=("/Users/adamcribbs/Documents/COMBATOxford/test_db/data/Pool_to_channel.tsv")
cursor.execute("create virtual table tmp46 using tsv("+sysdirs+")")


for ctime, channel in cursor.execute("select Pool, Channel from tmp45;"):
    print(ctime, channel)

    
