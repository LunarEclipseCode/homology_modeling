import os
import urllib.request
from urllib.error import HTTPError
from concurrent.futures import ThreadPoolExecutor, as_completed
import signal

alphabet = 'a'
start = '0'
end = 'z'

download_location = "D:\\Protein Database"

names = []

def main():
    name_creator(alphabet)
    file_checker()
    signal.signal(signal.SIGINT, handle_sigint)
    download_parallel()

def name_creator(letter):
    j = ord(start)
    while j <= ord(end):
        if j == ord(':'):
            j += 39
        for i in range(1, 10):
            k = ord('0')
            while k <= ord('z'):
                if k == ord(':'):
                    k += 39
                k += 1
                name = f"{i}{letter}{chr(j)}{chr(k)}"
                names.append(name)
        j += 1


def download_single(name):
    url = f"https://files.wwpdb.org/pub/pdb/data/structures/divided/pdb/{name[1:3]}/pdb{name}.ent.gz"
    destination = os.path.join(download_location, name[1], name[1:3], f"pdb{name}.ent.gz")
    # Create the directory if it doesn't exist
    os.makedirs(os.path.dirname(destination), exist_ok=True)
    try:
        urllib.request.urlretrieve(url, destination)
        return f"File {name} has been downloaded successfully."
    except HTTPError:
        return None
        
def download_parallel():
    with ThreadPoolExecutor() as executor:
        futures = {executor.submit(download_single, name): name for name in names}
        try:
            for future in as_completed(futures):
                name = futures[future]
                try:
                    result = future.result()
                    if result != None:
                        print(result)
                except Exception as e:
                    print(f"Error downloading file {name}: {str(e)}")
        
        except KeyboardInterrupt:
            print("Received KeyboardInterrupt. Shutting down gracefully.")
            
            # Cancel any remaining tasks
            for future in futures:
                future.cancel()
                
            executor.shutdown(wait=True)

def handle_sigint(signum, frame):
    print("Received SIGINT. Stopping execution.")
    os._exit(1)

def file_checker():
    c = 0
    while c < len(names):
        temp_file = os.path.join(download_location, names[c][1], names[c][1:3], f"pdb{names[c]}.ent.gz")
        if os.path.exists(temp_file):
            print(f"File {names[c]} already exists.")
            names.pop(c)
            c -= 1
        else:
            c += 1

if __name__ == "__main__":
    main()