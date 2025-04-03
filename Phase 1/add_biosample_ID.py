import csv
import concurrent.futures
import time
import random
import requests
import threading
import xml.etree.ElementTree as ET  # Import XML parsing library
from requests.exceptions import RequestException
import os # Import os to read environment variables

# Global rate limiter
MAX_REQUESTS_PER_MINUTE = 15  # Adjust as needed, erring on the side of caution
request_semaphore = threading.Semaphore(MAX_REQUESTS_PER_MINUTE)
last_reset_time = time.time()
request_count = 0

def fetch_assembly_details(bioproject_id, max_retries=5):
    global request_count, last_reset_time

    if not bioproject_id:
        return []

    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    api_key = os.environ.get("NCBI_API_KEY")  # Get API key from environment variable
    if not api_key:
        raise ValueError("NCBI_API_KEY environment variable not set.")

    for attempt in range(max_retries):
        try:
            # Global rate limiting
            request_semaphore.acquire()

            # Check if rate limit reset is needed
            current_time = time.time()
            if current_time - last_reset_time >= 60:
                with request_semaphore:  # Acquire all permits to reset
                    request_semaphore._value = MAX_REQUESTS_PER_MINUTE # Reset semaphore
                    request_count = 0
                    last_reset_time = current_time

            request_count += 1


            # Exponential backoff with a longer initial delay
            time.sleep((2 ** attempt) + random.random() + 5)  # Increased initial delay

            # Esearch request
            esearch_params = {
                "db": "assembly",
                "term": f"{bioproject_id}[BioProject]",
                "usehistory": "y",
                "retmax": 0,
                "api_key": api_key
            }
            esearch_response = requests.get(f"{base_url}esearch.fcgi", params=esearch_params)
            esearch_response.raise_for_status()

            # Extract WebEnv and query_key
            esearch_data = esearch_response.text
            web_env = esearch_data.split("<WebEnv>")[1].split("</WebEnv>")[0]
            query_key = esearch_data.split("<QueryKey>")[1].split("</QueryKey>")[0]

            # Efetch request
            efetch_params = {
                "db": "assembly",
                "query_key": query_key,
                "WebEnv": web_env,
                "rettype": "docsum",
                "retmode": "xml",
                "api_key": api_key
            }
            efetch_response = requests.get(f"{base_url}efetch.fcgi", params=efetch_params)
            efetch_response.raise_for_status()

            # Parse XML response
            try:
                root = ET.fromstring(efetch_response.content)  # Parse the XML content
                biosample_accns = [element.text for element in root.findall(".//BioSampleAccn")]
            except ET.ParseError as e:
                print(f"XML Parse Error: {e}")
                return []

            return biosample_accns

        except RequestException as e:
            print(f"Error on attempt {attempt + 1}: {str(e)}")
            if e.response is not None and e.response.status_code == 429:
                retry_after = e.response.headers.get('Retry-After')
                if retry_after:
                    print(f"Received Retry-After: {retry_after} seconds")
                    time.sleep(int(retry_after) + random.random())
                else:
                     time.sleep((2 ** attempt) + random.random() + 10)  # Longer delay on 429
            elif attempt == max_retries - 1:
                print(f"Failed to fetch data for BioProject {bioproject_id} after {max_retries} attempts")
                return []
        finally:
            # Correct semaphore release: Only release if it was acquired
            if request_semaphore._value <= MAX_REQUESTS_PER_MINUTE:
                request_semaphore.release()

# Function to process each row in parallel
def process_row(row):
    bioproject_id = row["BioProject_ID"]
    bio_sample_accns = fetch_assembly_details(bioproject_id)
    row["BioSampleAccn"] = bio_sample_accns
    return row, len(bio_sample_accns)  # Return row and count of BioSampleAccn values

def process_csv(input_csv, output_csv):
    with open(input_csv, 'r') as infile:
        reader = csv.DictReader(infile)
        rows = list(reader)

    max_biosample_count = 0
    all_biosample_data = []

    # Use ThreadPoolExecutor to speed up BioSampleAccn retrieval
    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as executor: # Reduced max_workers to 1
        futures = [executor.submit(process_row, row) for row in rows]
        for future in concurrent.futures.as_completed(futures):
            try:
                row, count = future.result()
                max_biosample_count = max(max_biosample_count, count)
                all_biosample_data.append(row)
            except Exception as e:
                print(f"Error processing row: {e}") # Catch potential errors from process_row

    # Generate column names dynamically
    fieldnames = reader.fieldnames + [f"BioSampleAccn_{i+1}" for i in range(max_biosample_count)]

    with open(output_csv, 'w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in all_biosample_data:
            biosample_accn_values = row.pop("BioSampleAccn", [])
            for i in range(max_biosample_count):
                row[f"BioSampleAccn_{i+1}"] = biosample_accn_values[i] if i < len(biosample_accn_values) else ""

            writer.writerow(row)

if __name__ == "__main__":
    input_csv = "/path/to/csv/with/bioprojectIDs"
    output_csv = "/path/to/output/csv"

    process_csv(input_csv, output_csv)

