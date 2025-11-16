# gw_agn_watcher/db.py
import requests
import psycopg2

def get_alerce_connection():
    """
    Fetch ALeRCE DB credentials and return a psycopg2 connection object.
    """
    url = "https://raw.githubusercontent.com/alercebroker/usecases/master/alercereaduser_v4.json"
    params = requests.get(url).json()['params']

    conn = psycopg2.connect(
        dbname=params['dbname'],
        user=params['user'],
        host=params['host'],
        password=params['password']
    )
    return conn
