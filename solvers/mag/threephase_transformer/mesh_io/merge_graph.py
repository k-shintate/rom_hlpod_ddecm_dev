from typing import Dict, List
import os
import argparse

def read_input_file(filename: str) -> Dict[int, List[int]]:
    """
    Reads the input file and returns a dictionary where:
      - The key is the first integer (index) on each line.
      - The value is a list of the remaining integers, excluding the second column.
    The first line of the file is skipped.
    """
    data: Dict[int, List[int]] = {}
    try:
        with open(filename, 'r') as file:
            next(file)  # Skip header line
            for line in file:
                numbers = list(map(int, line.split()))
                if len(numbers) < 2:
                    continue  # Skip line if not enough data
                index = numbers[0]      # Use the first number as the key
                values = numbers[2:]    # Skip the second column and take the rest as values
                if index in data:
                    data[index].extend(values)
                else:
                    data[index] = values
    except FileNotFoundError:
        print(f"Error: The file '{filename}' was not found.")
        print(f"Current working directory: {os.getcwd()}")
        raise
    return data

def write_output_file(data: Dict[int, List[int]], output_filename: str) -> None:
    """
    Writes the merged data to the output file.
    The format is:
      - First line: total number of entries.
      - Subsequent lines: each contains the index, the number of data points,
        and the data points separated by spaces.
    """
    with open(output_filename, 'w') as file:
        file.write(f"{len(data)}\n")
        # 安定した出力順にしたいなら sort
        for index in sorted(data.keys()):
            values = data[index]
            file.write(f"{index} {len(values)} {' '.join(map(str, values))}\n")

def merge_data(data1: Dict[int, List[int]], data2: Dict[int, List[int]]) -> Dict[int, List[int]]:
    """
    Merges two dictionaries containing lists of integers.
    If the same index exists in both, their lists are concatenated.
    """
    merged_data = data1.copy()
    for index, values in data2.items():
        if index in merged_data:
            merged_data[index] += values
        else:
            merged_data[index] = values
    return merged_data

def main(base_dir: str, elem_file: str, nedelec_file: str, output_file: str) -> None:
    input_filename1 = os.path.join(base_dir, elem_file)
    input_filename2 = os.path.join(base_dir, nedelec_file)
    output_filename = os.path.join(base_dir, output_file)

    try:
        data1 = read_input_file(input_filename1)
    except FileNotFoundError:
        print(f"Error reading file: {input_filename1}")
        return

    try:
        data2 = read_input_file(input_filename2)
    except FileNotFoundError:
        print(f"Error reading file: {input_filename2}")
        return

    merged_data = merge_data(data1, data2)
    write_output_file(merged_data, output_filename)
    print(f"Files merged successfully: {output_filename}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Merge graph_elem.dat and graph_nedelec_elem.dat into graph.dat in the given directory."
    )
    parser.add_argument(
        "base_dir",
        help="Base directory, e.g. result_mag/40-4-8"
    )
    parser.add_argument("--elem", default="graph_elem.dat", help="Input file name for elem (default: graph_elem.dat)")
    parser.add_argument("--nedelec", default="graph_nedelec_elem.dat", help="Input file name for nedelec (default: graph_nedelec_elem.dat)")
    parser.add_argument("--out", default="graph.dat", help="Output file name (default: graph.dat)")

    args = parser.parse_args()
    main(args.base_dir, args.elem, args.nedelec, args.out)
