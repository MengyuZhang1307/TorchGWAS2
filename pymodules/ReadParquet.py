import duckdb
import os
import time

def parquet_to_text_duckdb(
    input_file: str,
    output_file: str = None,
    num_threads: int = 16,
    memory_limit: str = "200GB",
):
    """
    Convert a Parquet file to a text file (CSV/TSV) using DuckDB efficiently.

    Parameters
    ----------
    input_file : str
        Path to the Parquet file.
    output_file : str, optional
        Path to save the text file. If None, adds '_converted.txt' to input name.
    num_threads : int, optional
        Number of CPU threads to use (default=8).
    memory_limit : str, optional
        Maximum memory allowed for DuckDB (default='200GB').
    """

    if output_file is None:
        output_file = os.path.splitext(input_file)[0] + "_converted.txt"

    start_time = time.time()

    # === Connect to DuckDB ===
    con = duckdb.connect(database=":memory:")
    con.execute(f"PRAGMA threads={num_threads};")
    con.execute(f"SET memory_limit='{memory_limit}';")
    con.execute("PRAGMA temp_directory='/tmp';")

    # === Build and run COPY query ===
    query = f"""
    COPY (
        SELECT * FROM '{input_file}'
    )
    TO '{output_file}'
    (DELIMITER '\t', HEADER, FORMAT 'csv');
    """

    print(f"Converting: {input_file}")
    print(f"Using {num_threads} threads")
    print(f"Memory limit: {memory_limit}")

    try:
        con.execute(query)
        print(f"\n Conversion complete!")
        print(f"Output saved to: {output_file}")
    except Exception as e:
        print(f"Error during conversion: {e}")
    finally:
        con.close()

    end_time = time.time()
    print(f"Time elapsed: {end_time - start_time:.2f} seconds")

