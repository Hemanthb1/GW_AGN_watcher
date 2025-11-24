# gw_agn_watcher/db.py
import requests
import psycopg2

import requests
import psycopg2

def get_alerce_connection():
    """
    Try connecting to ALeRCE DB using remote credentials.
    If that fails (fetch or connect), fall back to local parameters.
    """
    url = "https://raw.githubusercontent.com/alercebroker/usecases/master/alercereaduser_v4.json"

    # Define fallback parameters
    fallback_params = {
        "dbname": "alerce_local",
        "user": "alerceuser",
        "host": "localhost",
        "password": "your_fallback_password"
    }

    try:
        # --- Attempt remote fetch + connection ---
        params = requests.get(url, timeout=10).json()["params"]
        conn = psycopg2.connect(
            dbname=params["dbname"],
            user=params["user"],
            host=params["host"],
            password=params["password"]
        )
        print("✅ Connected to ALeRCE remote database.")
    except Exception as e:
        print(f"⚠️ Remote connection failed: {e}")
        print("🔁 Falling back to local database parameters...")
        conn = psycopg2.connect(
            dbname=fallback_params["dbname"],
            user=fallback_params["user"],
            host=fallback_params["host"],
            password=fallback_params["password"]
        )
        print("✅ Connected to local fallback database.")

    return conn
