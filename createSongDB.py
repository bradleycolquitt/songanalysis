#! /usr/bin/env python

import sys
import argparse
import MySQLdb as mdb


create_bird_table_sql = """CREATE TABLE IF NOT EXISTS birds\
                           (\
                               bird_id INT PRIMARY KEY,\
                               tags VARCHAR(8),\
                               manipulation VARCHAR(100),\
                               time_of_death DATETIME,\
                               song_data_dir VARCHAR(100),\
                               FOREIGN KEY id_birds_fk (bird_id)\
                                   REFERENCES sample_db.bird(idbirds)\
                           )
                        """
create_song_table_sql = '''CREATE TABLE IF NOT EXISTS songs\
                           (\
                               song_id INT AUTO_INCREMENT PRIMARY KEY,\
                               bird_id INT,\
                               song_datetime DATETIME,\
                               song_start FLOAT(3,2),\
                               song_end FLOAT(3,2),\
                               song_number SMALLINT,\
                               FOREIGN KEY bird_song_fk (bird_id)\
                                   REFERENCES birds(bird_id)\
                           )
                        '''

create_motif_table_sql = '''CREATE TABLE IF NOT EXISTS motifs\
                            (\
                                motif_id INT AUTO_INCREMENT PRIMARY KEY,\
                                song_id INT,\
                                motif_number TINYINT UNSIGNED,\
                                motif_start FLOAT(3,2),\
                                motif_end FLOAT(3,2),\
                                FOREIGN KEY song_motif_fk (song_id)\
                                    REFERENCES songs(song_id)\
                            )
                         '''

create_syllable_table_sql = '''CREATE TABLE IF NOT EXISTS syllables\
                               (\
                                   syllable_id INT AUTO_INCREMENT PRIMARY KEY,\
                                   motif_id INT,\
                                   syllable_number TINYINT UNSIGNED,\
                                   syllable_start FLOAT(3,2),\
                                   syllable_end FLOAT(3,2),\
                                   syllable_weinent FLOAT(5,5),\
                                   FOREIGN KEY motif_syllable_fk (motif_id)\
                                      REFERENCES motifs(motif_id)\
                               )
                            '''

create_gap_table_sql = '''CREATE TABLE IF NOT EXISTS gaps\
                          (\
                              gap_id INT AUTO_INCREMENT PRIMARY KEY,\
                              motif_id INT,\
                              gap_number TINYINT UNSIGNED,\
                              gap_duration SMALLINT UNSIGNED,\
                              FOREIGN KEY motif_syllable_fk (motif_id)\
                                  REFERENCES motifs(motif_id)\
                          )
                        '''

def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument(dest='db', help='Database name')
    parser.add_argument('-d', dest='drop', action='store_true', help='Drop tables if they exist')
    args = parser.parse_args()

    ## Create db
    conn = None
    cur = None
    try:
        conn = mdb.connect('localhost', 'brad', 'Eu23ler1')
        cur = conn.cursor()
        cur.execute("CREATE DATABASE IF NOT EXISTS {0}".format(args.db))
        cur.execute("USE {0}".format(args.db))
    except mdb.Error, e:
        print "MySQLdb error %d: %s " % (e.args[0], e.args[1])
        sys.exit()

    ## Create tables
    try:
        if args.drop:
            for x in ['gaps', 'syllables', 'motifs', 'songs', 'birds']:
                cur.execute("""DROP TABLES {0}""".format(x))
        cur.execute(create_bird_table_sql)
        cur.execute(create_song_table_sql)
        cur.execute(create_motif_table_sql)
        cur.execute(create_syllable_table_sql)
        cur.execute(create_gap_table_sql)
    except mdb.OperationalError, e:
        print "MySQLdb operational error %d: %s " % (e.args[0], e.args[1])

    conn.close()

if __name__ == "__main__":
    main(sys.argv)
