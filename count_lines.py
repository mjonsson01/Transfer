import os

def count_lines(directory):
    total_lines = 0
    file_count = 0

    # Extensions to include
    target_extensions = ('.h', '.cpp')

    print(f"{'File Name':<50} | {'Lines':<10}")
    print("-" * 65)

    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(target_extensions):
                file_path = os.path.join(root, file)
                try:
                    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                        lines = sum(1 for line in f)
                        total_lines += lines
                        file_count += 1
                        print(f"{file:<50} | {lines:<10}")
                except Exception as e:
                    print(f"Could not read {file}: {e}")

    print("-" * 65)
    print(f"Total Files: {file_count}")
    print(f"Total Lines of Code: {total_lines}")

if __name__ == "__main__":
    # Change '.' to your specific path if needed
    path_to_search = input("Enter the directory path (or press Enter for current): ") or "."
    count_lines(path_to_search)