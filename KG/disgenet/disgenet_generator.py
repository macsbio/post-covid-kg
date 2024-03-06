#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 11:29:24 2024

@author: alejandroadriaquelozano
"""
import sqlite3
import csv
import os 
from pathlib import Path
path = Path(__file__).resolve().parent
os.chdir(path)

def export_sqlite_to_csv(db_file, table_name, csv_file):
    # Connect to the SQLite database
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    # Execute a query to fetch all records from the specified table
    cursor.execute(f"SELECT * FROM {table_name}")

    # Fetch all rows from the cursor
    rows = cursor.fetchall()

    # Fetch the column names
    column_names = [description[0] for description in cursor.description]

    # Write the data to a CSV file
    with open(csv_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        # Write the column names as the header
        csv_writer.writerow(column_names)
        # Write the rows
        csv_writer.writerows(rows)

    # Close the database connection
    conn.close()
    
def get_table_names(db_file):
    # Connect to the SQLite database
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()

    # Query to fetch table names
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")

    # Fetch all table names from the cursor
    table_names = cursor.fetchall()

    # Close the database connection
    conn.close()

    return [name[0] for name in table_names]



# Example usage
db_file = 'disgenet_2020.db'  # SQLite database file
table_name = 'geneDiseaseNetwork'  # Table name within the database
csv_file = 'gene-disease_disgenet.csv'  # Output CSV file

table_names = get_table_names(db_file)
export_sqlite_to_csv(db_file, table_name, csv_file)
