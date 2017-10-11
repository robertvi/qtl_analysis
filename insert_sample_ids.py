#!/usr/bin/python

import MySQLdb

#main database
db = MySQLdb.connect(host="mongo",   #192.168.1.100
                     user="vicker",
                     passwd=open("/home/vicker/passwords/mysql_mongo_vicker").read(),
                     db="strawberry_samples")

cc = db.cursor()

for line in open('all_sample_ids'):
    sample_id = int(line.strip())
    cc.execute('insert into tmp_sample_list (sample_id) values (%d)'%sample_id)
    
db.commit()
db.close()
