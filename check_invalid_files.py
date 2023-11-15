import os

# checks for any irrelevant files in the protein database
def check_folder(directory):
    for root, dirs, files in os.walk(directory):
        for file in files:
            if not file.endswith('.gz'):
                # Print the directory if there is a file other than .gz
                print(f"Directory with non-.gz file: {root}")
                break  # Stop checking other files in the same directory

root_folder = "G:\Protein Database"
check_folder(root_folder)