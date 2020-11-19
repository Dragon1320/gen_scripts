import sqlite3

def dict_factory(cursor, row):
  d = {}
  for idx, col in enumerate(cursor.description):
    d[col[0]] = row[idx]
  return d

def get_db_cursor(db_path):
  conn = sqlite3.connect(db_path)
  conn.row_factory = dict_factory
  
  return conn.cursor()
